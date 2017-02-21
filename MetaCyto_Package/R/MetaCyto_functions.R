# Define fcsFrame functions----------------

#' label each clusters as "+" or "-" or neutral for each marker
#'
#' A function that label each cluster as "+" or "-" or neutral for each marker
#' @param fcsFrame A flowFrame object.
#' @param clusterList a list, each element should be a vector containing the IDs of all cells that are belong to a cluster
#' @param excludeClusterParameters A vector specifiying the name of markers not to be used for labeling.
#' @param minPercent a number between 0 and 0.5. Used to specify the minimum percent of cells in positive and negative region after bisection. Keep it small to avoid bisecting uni-mode distributions.
#' @param labelQuantile a number between 0 and 0.5. Used to specify the minimum percent of a cluster required to be larger or smaller than the cutoff value for labeling.
#' @return Returns a list with two components: 1) clusterLabel, contains a vector of lables, each correspond to a cluster in clusterList. 2) cutoff, contains a vector of cutoff values used to bisect each marker.
#' @export
labelCluster= function(fcsFrame,
                       clusterList,
                       excludeClusterParameters=c("TIME"),
                       minPercent=0.05,
                       labelQuantile=0.95,
                       cutoff=NULL){

  #prepare exclude parameters
  excludeClusterParameters=toupper(excludeClusterParameters)
  # get marker names
  antibodies=markerFinder(fcsFrame)

  # Get expression matrix
  expr=flowCore::exprs(fcsFrame);
  colnames(expr)=antibodies

  # subset on columns
  w = (!antibodies%in%excludeClusterParameters) #& (!channels%in%excludeClusterParameters)
  antibodies=antibodies[w]
  expr=expr[,w,drop=F]

  #find cutoff of each parameter
  if(is.null(cutoff)){
    cutoff=apply(X=expr,MARGIN=2,FUN=findCutoff)
  }

  #label each cluster
  CL_label=rep(NA,length(clusterList))
  no_label=sapply(1:length(antibodies),function(i){
    x1=quantile(expr[,i],minPercent)
    x2= quantile(expr[,i],(1-minPercent))
    if(cutoff[i]<x1|cutoff[i]>x2){return(T)}else{return(F)}
  })
  for(i in 1:length(clusterList)){
    l=c()
    for(j in 1:length(antibodies)){
      if(no_label[j]==T){next}
      x=expr[clusterList[[i]],j]
      if(quantile(x,labelQuantile)<cutoff[j]){
        l=c(l,paste0(antibodies[j],"-"))
      }else if(quantile(x,(1-labelQuantile))>cutoff[j]){
        l=c(l,paste0(antibodies[j],"+"))
      }
    }
    l=paste(l,collapse="|")
    CL_label[i]=l
  }

  # merge clusters with the same labels
  CL_label=unique(CL_label)

  # find cluster with 1 sign difference, merge to create new cluster
  while(length(CL_label)>1){
    label_pair=combn(CL_label,2)
    new_label=apply(label_pair,2,labelCombiner)
    w=which(new_label%in%CL_label)
    new_label[w]=NA
    if(sum(!is.na(new_label))==0){break}
    w=which(!is.na(new_label))
    CL_label=c(CL_label,new_label[w])
  }

  CL_label=labelUnifier(CL_label)
  CL_label=unique(CL_label)
  Result=list("clusterLabel"=CL_label,"cutoff"=cutoff)
  return(Result)
}

#' search for clusters using pre-defined labels
#'
#' A function that search for clusters using pre-defined labels
#' @param fcsFrame A flowFrame object.
#' @param clusterLabel a vector of labels, such as "CD3+|CD4+|CD8-"
#' @param cutoff a vector of cutoff values to bisect the distribution of each
#'   marker. The names of the vector should be the same as the marker names. If
#'   NULL, the cutoff value will be determined automatically.
#' @param rmNULL True or False. Used to specifiy if a cluster with 0 cells
#'   should be returned or not.
#' @return Returns a list with two components: 1) clusterList, contains a list.
#'   Each element of the list is a vector containing the ID of all cells in a
#'   cluster. The names corresponds to the labels specified in clusterLabel. 2)
#'   cutoff, contains a vector of cutoff values used to bisect each marker.
#' @export
searchCluster=function(fcsFrame,
                       clusterLabel,
                       cutoff=NULL,
                       rmNull=T){

  #prepare exclude parameters
  clusterLabel=toupper(clusterLabel)
  antibodies=markerFinder(fcsFrame)

    # Get expression matrix
  expr=flowCore::exprs(fcsFrame);
  colnames(expr)=antibodies
  AB=unlist(strsplit(clusterLabel,split="&|\\|"))
  AB=sapply(AB,function(x){gsub("\\+$|-$|\\^NE$|\\^LO$|\\^HI$","",x)})
  w=which(antibodies%in%AB)
  expr=expr[,w,drop=F]

  #find labels whose markers are included in the study
  CL_label=clusterLabel
  CL_label2=strsplit(clusterLabel,split="&|\\|")
  in_study=sapply(CL_label2,function(x){
    Ab_in_label=gsub("\\+$|-$|\\^NE$|\\^LO$|\\^HI$","",x)
    if(length(setdiff(Ab_in_label,antibodies))>0){return(F)}else{return(T)}
  })
  CL_label=CL_label[in_study]
  CL_label2=CL_label2[in_study]
  if(length(CL_label)<1){return(list())}

  #find cutoff of each parameter
  if(is.null(cutoff)){
    cutoff=apply(X=expr,MARGIN=2,FUN=findCutoff)
  }

  #find each cluster
  CL=list()
  for(i in 1:length(CL_label2)){
    rows=sapply(CL_label2[[i]],function(x){
      m=gsub("\\+$|-$|\\^NE$|\\^LO$|\\^HI$","",x)
      if(grepl("\\+$",x)){w=which(expr[,m]>cutoff[m])}
      if(grepl("-$",x)){w=which(expr[,m]<cutoff[m])}
      if(grepl("\\^NE$",x)){w=which(expr[,m]<triS[2,m])}
      if(grepl("\\^LO$",x)){w=which(expr[,m]<triS[1,m])&expr[,m]>triS[2,m]}
      if(grepl("\\^HI$",x)){w=which(expr[,m]>triS[1,m])}
      return(w)
    })
    if(length(CL_label2[[i]])>1){rows=Reduce(intersect,rows)}
    CL[[i]]=rows
  }
  if(rmNull==T){
    CL_count=sapply(CL,length)
    CL=CL[CL_count>1]
    CL_label=CL_label[CL_count>1]
  }

  if(length(CL)<1){return(list())}

  names(CL)=CL_label
  Result=list("clusterList"=CL,"cutoff"=cutoff)
  return(Result)
}

#' preprocess fcs files from a single experiment
#'
#' A function that preprocess fcs files from a single experiment.
#' @param fcsFiles A vector specifying the location of all fcs files.
#' @param assay Either "FCM" or "CyTOF" to indicate the type of cytometry data.
#' @param b a positive number used to specify the arcsinh transformtion. f(x) = asinh (b*x) where x is the original value and f(x) is the value after transformation. The suggested value is 1/150 for flow cytometry (FCM) data and 1/8 for CyTOF data.
#' @param fileSampleSize An integer specifying the number of events sampled from each fcs file. If NULL, all the events will be pre-processed and wrote out to the new fcs files.
#' @param excludeTransformParameters A vector specifiying the name of parameters not to be transformed (left at linear scale).
#' @return Returns a flowFrame object containing the preprocessed cytometry data. Cells from different fcs files are combined into one flow frame. A new paramter, sample_id, is introduced to indicate the origin of each cell.
#' @export
preprocessing=function(fcsFiles,
                       assay=c("FCM", "CyTOF"),
                       b=1/200,
                       fileSampleSize=5000,
                       excludeTransformParameters=c("FSC-A","FSC-W","FSC-H","Time","Cell_length")){

  fcs_param=NULL
  fcs_names= gsub(".*/","",fcsFiles)
  excludeTransformParameters=paste(excludeTransformParameters,collapse="|")


  cat("Preprocessing \n")

  # 1) identify sample :
  fcsFiles=as.character(fcsFiles)

  # 2) read the fcs files
  fcs = flowCore::read.flowSet(fcsFiles,transformation="linearize",alter.names=FALSE,truncate_max_range=F)
  if(!is.null(fileSampleSize)){
    fcs = flowCore::fsApply(fcs, function(f){
      L=nrow(f@exprs)
      if(L>fileSampleSize){
        f@exprs=f@exprs[sample(1:L,fileSampleSize),]
      }
      return(f)
    })
  }

  # 3) compensation and transformation
  if (assay == "FCM") {
    # retreiving fluorescent biomarkers
    w=which(!grepl(excludeTransformParameters,flowCore::colnames(fcs),ignore.case = T))
    biomarker_vector = flowCore::colnames(fcs)[w]
    # identify the fcs file format
    version = flowCore::fsApply(fcs, function(frame) {
      return(flowCore::keyword(frame, "FCSversion")[[1]])
    })
    unique_version = as.numeric(unique(as.vector(version)))

    # compensation, transformation
    ## if the file version is 2.0, then code will log transform the data
    if (unique_version == 2) {
      trans = flowCore::logTransform()
      translist = flowCore::transformList(biomarker_vector,trans)
      fcs = flowCore::transform(fcs, translist)
    } else if (unique_version == 3) {

      ### function first checks if the spill matrix in the flowframe is an actual spill matrix or an identity matrix
      if (is.null(flowCore::keyword(fcs[[1]], "SPILL")[[1]]) == FALSE) {
        check = flowCore::fsApply(fcs, function(x) {
          result = isSymmetric(flowCore::keyword(x, "SPILL")[[1]])
        })
        ### if all flowframes have spill matrices then my function
        ### will apply them to each flowframe in the flowset for compensation
        if (unique(check)[[1]] == FALSE) {
          fcs = flowCore::fsApply(fcs, function(x) {
            new_frame = flowCore::compensate(x, flowCore::keyword(x, "SPILL")[[1]])
            return(new_frame)
          })
        }
      }

      trans = flowCore::arcsinhTransform(transformationId="defaultArcsinhTransform",a=0,b=b,c=0)
      translist = flowCore::transformList(biomarker_vector, trans)
      fcs = flowCore::transform(fcs, translist)
    }
  } else if (assay == "CyTOF") {
    w=which(!grepl(excludeTransformParameters,flowCore::colnames(fcs),ignore.case = T))
    biomarker_vector = flowCore::colnames(fcs)[w]

    trans = flowCore::arcsinhTransform(transformationId = "defaultArcsinhTransform", a = 0, b = b, c = 0)
    translist = flowCore::transformList(biomarker_vector, trans)
    fcs = flowCore::transform(fcs, translist)
  }
  fcs=set2Frame(fcs)
  return(fcs)
}

#' derive summary statistics for clusters
#'
#' A function that derives summary statistics for clusters.
#' @param fcsFrame A flow Frame object returned from preprocessing function. Must contain a parameter called "sample_id" that specifiy the origin of each cell using integer IDs.
#' @param clusterList a list, each element should be a vector containing the IDs of all cells that are belong to a cluster
#' @param fcsNames a vector of fcs file names. Each element corresponds to an interger ID in the "sample_id" parameter in fcsFrame.
#' @return Returns a data frame, with rows correspond to each fcs file, columns correspond to markers or fractions.
#' @export
clusterStats=function(fcsFrame,clusterList,fcsNames){
  antibodies=markerFinder(fcsFrame)
  expr=flowCore::exprs(fcsFrame)
  colnames(expr)=antibodies
  CL=clusterList
  CL_label=names(clusterList)
  sample_name=fcsNames
  expr_sample=expr[,"SAMPLE_ID"]
  expr=expr[,colnames(expr)!="SAMPLE_ID"]
  antibodies=antibodies[antibodies!="SAMPLE_ID"]
  SP_stat=NULL
  for(i in 1:length(CL)){
    for(j in unique(expr_sample)){
      fileSampleSize=sum(expr_sample==j)
      w=intersect( which(expr_sample==j), CL[[i]] )
      frac=length(w)/fileSampleSize
      if(length(w)>1){M=apply(expr[w,],2,median)}
      if(length(w)==1){M=expr[w,]}
      if(length(w)==0){M=rep(NA,ncol(expr))}
      t1=c(paste0("cluster",i),CL_label[i],sample_name[j],M,frac)
      SP_stat=rbind(SP_stat,t1)
    }
  }
  colnames(SP_stat)=c("cluster_id","label","fcs_names",antibodies,"fraction")
  rownames(SP_stat)=NULL
  SP_stat=data.frame(SP_stat,stringsAsFactors = F,check.names = F)
  SP_stat[,1:3]=lapply(SP_stat[,1:3],as.character)
  SP_stat[,-c(1:3)]=lapply(SP_stat[,-c(1:3)],as.numeric)

  return(SP_stat)
}

#' draw density plot for each cluster.
#'
#' A function that draw density plot for each cluster.
#' @param fcsFrame A flow Frame object returned from preprocessing function.
#' @param clusterList a list, each element should be a vector containing the IDs of all cells that are belong to a cluster
#' @param cutoff a vector of cutoff values to bisect the distribution of each
#'   marker. The names of the vector should be the same as the marker names.
#' @param markerToPlot a vector specifying markers included in the plot. If NULL, all markers will be plotted.
#' @return NULL. The plot will show up automatically.
#' @export
densityPlot=function(fcsFrame,clusterList,cutoff,markerToPlot=NULL){
  CL_label=toupper(names(clusterList))
  markerToPlot=toupper(markerToPlot)
  antibodies=markerFinder(fcsFrame)
  expr=flowCore::exprs(fcsFrame);
  colnames(expr)=antibodies
  if(length(markerToPlot)<1){markerToPlot=antibodies}
  #antibodies=toupper(names(cutoff))
  #expr=expr[,antibodies]
  CL=clusterList
  par( mfcol = c(length(CL), length(markerToPlot) ) )
  for(i in 1:length(markerToPlot)){
    x_all=expr[,markerToPlot[i]]
    for(j in 1:length(CL)){
      b=seq(min(x_all),max(x_all), ((max(x_all)-min(x_all))/100) )
      hist(x_all,col=rgb(0, 0, 0, 0.2),xlab=markerToPlot[i],breaks=b,freq=T,border=F,main=paste0(markerToPlot[i]," of cluster ",j," : \n", CL_label[j]))
      hist(expr[CL[[j]],markerToPlot[i]],add=T,breaks=b,col=rgb(1, 0, 0, 0.5),freq=T,border=F)
      if(markerToPlot[i]%in%names(cutoff)){abline(v=cutoff[markerToPlot[i]])}
    }
  }
}


# Define batch functions----------------

#' preprocessing fcs files from different studies in batch.
#'
#' It transform and compensate for the raw fcs files and write out the processed data to a new set of fcs files.
#' @param inputMeta A dataframe containing 2 columns: a column called "fcs_files" that contains the location (relative to the working directory) of each fcs file on the hard drive and a column called "study_id" that specify what study each fcs file belongs to.
#' @param assay Either "FCM" or "CyTOF" to indicate the type of cytometry data.
#' @param outpath a string indicating the directory the pre-processed fcs files will be wrote to.
#' @param b a positive number used to specify the arcsinh transformtion. f(x) = asinh (b*x) where x is the original value and f(x) is the value after transformation. The suggested value is 1/150 for flow cytometry (FCM) data and 1/8 for CyTOF data.
#' @param fileSampleSize An integer specifying the number of events sampled from each fcs file. If NULL, all the events will be pre-processed and wrote out to the new fcs files.
#' @param excludeTransformParameters A vector specifiying the name of parameters not to be transformed (left at linear scale).
#' @return Does not return anything. The output is wrote to the directory specified by the "outpath".
#' @details The function takes a data frame which specify the location of the fcs files and the panels the fcs files belong to. It transform the cytometry data using the arcsinh transformtion. For flow cytometry data, it compensate the data using the compensation matrix supplied in the fcs file. the preprocessed fcs files and a table called "processed_sample_summary.csv" will be wrote out to outpath as well.
#' @export

preprocessing.batch = function(inputMeta,
                               assay=c("FCM", "CyTOF"),
                               outpath,
                               b=1/150,
                               fileSampleSize=5000,
                               excludeTransformParameters=c("FSC-A","FSC-W","FSC-H","Time","Cell_length")){

  fcs_param=NULL
  fcs_names= gsub(".*/","",inputMeta$fcs_files)
  excludeTransformParameters=paste(excludeTransformParameters,collapse="|")
  dir.create(outpath,recursive=T)
  if(length(b)==1){b=rep(b,nrow(inputMeta))}
  if(length(assay)==1){assay=rep(assay,nrow(inputMeta))}
  for(std in unique(inputMeta$study_id)){

    cat("Study ID = ",std, " ")

    # 1) identify sample :
    fcs_files=subset(inputMeta$fcs_files,inputMeta$study_id==std)
    fcs=preprocessing(fcsFiles=fcs_files,
                      assay=subset(assay,inputMeta$study_id==std)[1],
                      b=subset(b,inputMeta$study_id==std)[1],
                      fileSampleSize=fileSampleSize,
                      excludeTransformParameters=excludeTransformParameters)

    # 5) outputting resuts
    flowCore::write.FCS(fcs,filename=paste0(outpath,'/',std,".fcs"))

    antibodies=markerFinder(fcs)
    antibodies=paste(antibodies,collapse="|")
    fcs_param=rbind(fcs_param,data.frame("fcs_files"=fcs_files,"study_id"=std,"antibodies"=antibodies))
  }
  fcs_param=cbind("fcs_names"=fcs_names,fcs_param)
  write.csv(fcs_param,paste0(outpath,"/processed_sample_summary.csv"),row.names=F)
  cat("Preprocess result stored in the folder:",outpath, "\n")
}#end of function


#' cluster the pre-processed fcs files from different studies in batch
#'
#' A function that clusters the pre-processed fcs files from different studies in batch
#' @param preprocessOutputFolder directory where the pre-processed results are stored.
#' @param excludeClusterParameters A vector specifiying the name of markers not to be used for clustering and labeling.
#' @param labelQuantile a number between 0 and 0.5. Used to specify the minimum percent of a cluster required to be larger or smaller than the cutoff value for labeling.
#' @param clusterFunction The name of unsuperviesed clustering function the user wish to use for clustering the cells. The default is "flowSOM.MC". The first argument of the function must take a flow frame, the second argument of the function must take a vector of excludeClusterParameters. The function must returns a list of clusters containing cell IDs. flowSOM.MC and flowHC are implemented in the package. For other methods, please make your own wrapper functions.
#' @param minPercent a number between 0 and 0.5. Used to specify the minimum percent of cells in positive and negative region after bisection. Keep it small to avoid bisecting uni-mode distributions.
#' @param ... pass arguments to labelCluster and clusterFunction
#' @return a vector of labels identified in the cytometry data.
#' @export
autoCluster.batch= function(preprocessOutputFolder,
                            excludeClusterParameters=c("TIME"),
                            labelQuantile=0.9,
                            clusterFunction=flowSOM.MC,
                            minPercent=0.05,...){
  #read the output from preprocessing
  inputMeta=read.csv(paste0(preprocessOutputFolder,'/processed_sample_summary.csv'),stringsAsFactors=F)
  #create output foler

  #prepare exclude parameters
  excludeClusterParameters=toupper(excludeClusterParameters)
  all_labels=NULL
  for(std in unique(inputMeta$study_id)){
    cat("Clustering , study ID = ",std, "\n")

    ##### 1) read sample files for each study############################################################
    fcs_files=paste0(preprocessOutputFolder,"/",std,".fcs")
    fcs=flowCore::read.FCS(fcs_files,truncate_max_range=F)

    # make sure the fcs file antibody names are the same as the preprocessed output
    antibodies=subset(inputMeta$antibodies,inputMeta$study_id==std)[1]
    antibodies=strsplit(antibodies,"\\|")[[1]]

    ##### 2) subset the cells in fcs ############################################################
    # Get expression matrix
    expr=flowCore::exprs(fcs);
    colnames(expr)=antibodies

    # subset on columns
    w=!antibodies%in%excludeClusterParameters
    antibodies=antibodies[w]
    if(length(antibodies)<2){next}
    expr=expr[,w,drop=F]
    expr_scale=scale(expr,center=F,scale=T)
    fcs=flowCore::flowFrame(expr_scale)

    CL=clusterFunction(fcs,excludeClusterParameters,...)
    CL_label=labelCluster(fcs,CL,excludeClusterParameters,
                          labelQuantile=labelQuantile,
                          minPercent=minPercent,...)
    all_labels=union(all_labels,CL_label$clusterLabel)
  }#end of each study
  return(all_labels)
}

#' search for clusters using pre-defined labels in cytometry data from different studies in batch
#'
#' A function that search for clusters using pre-defined labels in cytometry data from different studies in batch.
#' @param preprocessOutputFolder directory where the pre-processed results are stored.
#' @param outpath a string indicating the directory the results should be write to.
#' @param clusterLabel a vector containing labels, such as "CD3+|CD4+|CD8-"
#' @param ifPlot True or False. Used to specifiy if a the density plot for each cluster should be plotted
#' @return NULL
#' @details The function write out the summary statistics for each cluster. A seperate directory will be created for each study.
#' @export
searchCluster.batch=function(preprocessOutputFolder,
                             outpath="search_output",
                             clusterLabel,
                             ifPlot=T){

  #read the output from preprocessing
  inputMeta=read.csv(paste0(preprocessOutputFolder,'/processed_sample_summary.csv'),stringsAsFactors=F)
  #create output foler
  dir.create(file.path(outpath),recursive=T)

  #prepare exclude parameters
  clusterLabel=toupper(clusterLabel)
  cluster_summary=NULL
  for(std in unique(inputMeta$study_id)){
    cat("Searching , study ID = ",std, "\n")
    dir.create(paste(outpath,std,sep="/"))
    ##### 1) read sample files for each study############################################################
    fcs_files=paste0(preprocessOutputFolder,"/",std,".fcs")
    fcs=flowCore::read.FCS(fcs_files,truncate_max_range=F)

    # make sure the fcs file antibody names are the same as the preprocessed output
    antibodies=subset(inputMeta$antibodies,inputMeta$study_id==std)[1]
    antibodies=strsplit(antibodies,"\\|")[[1]]

    ##### 2) subset the cells in fcs ############################################################
    # Get expression matrix
    expr=flowCore::exprs(fcs);
    colnames(expr)=antibodies
    fcs=flowCore::flowFrame(expr)

    ##### 3) bisect each marker############################################################
    sC=searchCluster(fcsFrame=fcs,clusterLabel=clusterLabel)
    CL=sC$clusterList
    CL_label=names(CL)
    if(length(CL)<1){next}

    CL_stats=clusterStats(fcs,CL,inputMeta$fcs_names[inputMeta$study_id==std])
    CL_stats=cbind("study_id"=std,"fcs_files"=inputMeta$fcs_files[inputMeta$study_id==std],CL_stats)
    write.csv(CL_stats,paste(outpath,std,"cluster_stats_in_each_sample.csv",sep="/"),row.names=F)

    ##### 6) plot the density plot############################
    if(ifPlot==T){
      filename=paste(outpath,std,"density_plot.pdf",sep="/")
      height=3*length(CL);width=3*length(antibodies)
      pdf(filename,width=width,height=height)
      densityPlot(fcs,CL,cutoff=sC$cutoff)
      dev.off()
    }
  }#end of each study
}
# Define supporting functions---------

#' organize fcs files in a study from ImmPort into panels
#'
#' A function that organizes fcs files in a study from ImmPort into panels.
#' @param metaData a data frame. Must contain a column listing the names of fcs files included in the study.
#' @param studyFolder Path of directory containing all the files of a study from ImmPort.
#' @param fcsCol a string specifiying the name of the column in metaData that lists fcs files included in the study.
#' @param assay Either "FCM" or "CyTOF" to indicate the type of cytometry data.
#' @return A dataframe containing 2 columns: a column called "fcs_files" that contains the location (relative to the working directory) of each fcs file on the hard drive and a column called "study_id" that specify what study each fcs file belongs to.
#' @export
fcsInfoParser =function(metaData,
                        studyFolder,
                        fcsCol="ZBXFN",
                        assay=c("FCM", "CyTOF")){
  cat("Parsing the meta-data from ImmPort\n")
  metaData=unique(metaData[,fcsCol,drop=FALSE])

  if(assay=="FCM"){
    fcs_files=paste(studyFolder,"ResultFiles/Flow_cytometry_result",metaData[,fcsCol],sep="/")
  }else{
    fcs_files=paste(studyFolder,"ResultFiles/CyTOF_result",metaData[,fcsCol],sep="/")
  }
  fsc_markers=sapply(fcs_files,function(x){
    tryCatch({
      fcs=flowCore::read.FCS(x,truncate_max_range=F)
      panel=paste(flowCore::pData(flowCore::parameters(fcs))$name,
                  flowCore::pData(flowCore::parameters(fcs))$desc, sep="-")
      return(paste(panel,collapse="|"))
    },error=function(e){return("DoesNotExist")})
  })
  study=sapply(fsc_markers,function(x){
    if(x=="DoesNotExist" ){return("DoesNotExist")}
    w=which(unique(fsc_markers)==x);
    w=paste0(gsub("/","_",studyFolder),"_",assay,"-",w)
    return(w)
  })
  names(study)=NULL
  inputMeta=data.frame("fcs_files"=fcs_files,"study_id"=study)

  inputMeta=subset(inputMeta,!grepl("DoesNotExist",inputMeta$study_id))
  inputMeta=unique(inputMeta)
  inputMeta$fcs_files=as.character(inputMeta$fcs_files)
  return(inputMeta)
}

#' collect and combine data from multiple csv files of the same format
#'
#' A function that collect and combine data from multiple csv files of the same format.
#' @param files a vector containing the paths of csv files to be combined.
#' @param longform True or False. Used to specify if the table in each csv file should be converted into long form before combining.
#' @return A dataframe containing conbined information from multiple csv files.
#' @export
collectData=function(files,longform=T){
  all_data=NULL
  nms=NULL
  for(fn in files){
    #read output data from "autoCluster" function
    cluster_stat=read.csv(fn,stringsAsFactors=F,check.names=F)
    if(is.null(nms)){nms=colnames(cluster_stat)}
    w=which(nms=="fcs_names")+1
    if(longform==T){
      cluster_stat=tidyr::gather(cluster_stat,parameter_name, value ,w:ncol(cluster_stat))
    }else(colnames(cluster_stat)=nms)
    all_data=rbind(all_data,cluster_stat)
  }
  return(all_data)
}

#' collect sample information for fcs files in  a study from ImmPort.
#'
#' A function that collect sample information for fcs files in  a study from ImmPort..
#' @param metaData a data frame. Must contain a column listing the names of fcs files included in the study.
#' @param studyFolder Path of directory containing all the files of a study from ImmPort.
#' @param fcsCol a string specifiying the name of the column in metaData that lists fcs files included in the study.
#' @param assay Either "FCM" or "CyTOF" to indicate the type of cytometry data.
#' @param attrCol a vector of column names. Used to specify the information about each cytometry the user wish to include in the analysis.
#' @return A dataframe containing sample information.
#' @export
sampleInfoParser =function(metaData=selected_data,
                           studyFolder,
                           fcsCol="ZBXFN",
                           assay="FCM",
                           attrCol){
  if(assay=="FCM"){
    fcs_files=paste(studyFolder,"ResultFiles/Flow_cytometry_result",metaData[,fcsCol],sep="/")
  }else{
    fcs_files=paste(studyFolder,"ResultFiles/CyTOF_result",metaData[,fcsCol],sep="/")
  }
  inputMeta=data.frame("fcs_files"=fcs_files)
  inputMeta=cbind(inputMeta,metaData[,attrCol])
  inputMeta=unique(inputMeta)
  inputMeta$fcs_files=as.character(inputMeta$fcs_files)
  return(inputMeta)
}


#' summarize markers in panels.
#'
#' A function that summarize markers in panels.
#' @param panelInfo a data frame returned by collectData function. It should contain all the information outputted by the preprocessing.batch function.
#' @param folder the directory where the output should be written
#' @param cluster True or False. Used to indicate if the markers and panels should be clustered in the plot.
#' @param plotImage True or False. Used to indicate if a plot summarizing markers in  panels should be produced.
#' @param width Used to specify the width of the plot
#' @param height Used to specify the height of the plot
#' @return A dataframe describing what markers are in each panel.
#' @export
panelSummary=function(panelInfo,folder,cluster=T,plotImage=T,width=20,height=20){
  panelInfo=unique(panelInfo[,c("study_id","antibodies")])
  ab_list=NULL
  for(i in 1:nrow(panelInfo)){
    t1=strsplit(panelInfo$antibodies[i],split="\\|")[[1]]
    t1=t1[t1!="NA"&is.na(t1)==F]
    if(length(t1)>0){
      ab_list=rbind(ab_list, data.frame("study_id"=panelInfo$study_id[i],"antibodies"=t1,"value"=1))
    }
  }
  ab_list=unique(ab_list)
  #ab_list=unique(ab_list)
  ab_table=tidyr::spread(ab_list,key=antibodies,value=value,fill=0)
  rownames(ab_table)=ab_table$study_id;ab_table=ab_table[,-1]
  ab_table=t(data.matrix(ab_table))
  if(cluster==T){
    d=as.dist(1-cor(ab_table))
    c1=hclust(d)$order
    d=as.dist(1-cor(t(ab_table)))
    r1=hclust(d)$order
    ab_table=ab_table[r1,c1]
  }
  write.csv(ab_table,paste0(folder,"/panel_summary.csv"))

  if(plotImage==T){
    pdf(paste0(folder,"/panel_summary.pdf"),width=width,height=height)
    op <- par(mar = c(30,30,30,10))
    image(z = ab_table, col = c("white","red"), axes = F)
    t1=1/(nrow(ab_table)-1)
    axis(side = 3, labels = rownames(ab_table),las=2,cex.axis=3,
         at = seq(0, 1,by = t1))
    t1=1/(ncol(ab_table)-1)
    axis(side = 2, labels = colnames(ab_table),las=2,cex.axis=3,
         at = seq(0, 1,by = t1))
    box()
    par(op)
    dev.off()
  }
  return(ab_table)
}




#' Used to update marker names
#'
#' A function that update marker names in the files output by the preprocessing.batch function.
#' @param oldNames a vector of marker names you wish to change
#' @param newNames a vector of marker names you wish each oldNames to be changed to.
#' @param files a list of "processed_sample_summary.csv" files in which the name change will ocurre.
#' @return Null
#' @export
nameUpdator=function(oldNames,newNames,files){
  oldNames=as.vector(oldNames)
  newNames=as.vector(newNames)
  for(fn in files){
    TB=read.csv(fn,stringsAsFactors=F,check.names=F)
    for(i in 1:length(newNames)){
      TB$antibodies=sapply(TB$antibodies,function(x){
        gsub(oldNames[i],newNames[i],x,fixed=T)
      })
    }
    write.csv(TB,fn,row.names=F)
  }
}

#' combine cells in a flow set into a flow frame.
#'
#' A function that combine cells in a flow set into a flow frame.
#' @param flowSet a flow set object
#' @return Returns a flowFrame object. All cells from flow set are combined into one flow frame. A new paramter, sample_id, is introduced to indicate the origin of each cell.
#' @export
set2Frame=function(flowSet){
  expr=flowCore::fsApply(flowSet,function(x){
    v=flowCore::exprs(x);
    return(v)
  })
  eventN=flowCore::fsApply(flowSet,function(x){
    v=flowCore::exprs(x);
    n=nrow(v)
    return(n)
  })
  Label=flowCore::fsApply(flowSet,function(x){
    v=flowCore::pData(flowCore::parameters(x));
    return(v)
  })
  #annotate expr colnames
  channels=Label[[1]]$name
  antibodies=Label[[1]]$desc
  antibodies=toupper(antibodies)
  channels=toupper(channels)
  antibodies=sapply(1:length(antibodies), function(i){
    if(is.na(antibodies[i])|antibodies[i]=="NA"){return(channels[i])}else{antibodies[i]}
  })
  colnames(expr)=antibodies
  #add sample_id
  sample_id=lapply(1:length(eventN),function(i){rep(i,eventN[i])})
  sample_id=unlist(sample_id)
  expr=cbind(expr,"sample_id"=sample_id)
  #return flowFrame
  fFrame=flowCore::flowFrame(expr)
  return(fFrame)
}

#' filter cluster labels
#'
#' A function that filter cluster labels.
#' @param labels a vector containing labels for cell clusters
#' @param minPlus an integer, used to specify the minmum number of "+" a label should contain.
#' @param minMarker an integer, used to specify the minmum number of markers a label should contain.
#' @param maxMarker an integer, used to specify the max number of markers a label should contain.
#' @return returns a vector of labels that pass through the filter.
#' @export
filterLabels=function(labels,minPlus,minMarker,maxMarker){
  plus_N=sapply(labels,function(x){
    lengths(regmatches(x, gregexpr("\\+", x)))
  })
  marker_N=sapply(labels,function(x){
    lengths(regmatches(x, gregexpr("\\|", x)))+1
  })
  labels=labels[plus_N>=minPlus&marker_N<=maxMarker&marker_N>=minMarker]
  return(labels)
}

#' find markers in a flow frame object
#'
#' A function that finds markers in a flow frame object.
#' @param fcsFrame a flow frame object.
#' @return returns a vector of markers.
#' @details If the antibody name is available, the antibody name will be returned, otherwise the channel name will be returned.
#' @export
markerFinder=function(fcsFrame){
  Label=flowCore::pData(flowCore::parameters(fcsFrame));
  channels=Label$name
  antibodies=Label$desc
  antibodies=toupper(antibodies)
  channels=toupper(channels)
  antibodies=sapply(1:length(antibodies), function(i){
    if(is.na(antibodies[i])|antibodies[i]=="NA"){return(channels[i])}else{antibodies[i]}
  })
  return(antibodies)
}

# Define meta-analysis functions----------

#' perform meta-analysis
#'
#' A function that performs meta-analysis
#' @param value a string to specify the column name of the dependent variable (y)
#' @param variableOfInterst a string to specify the column name of the independent variable of interest (x1)
#' @param otherVariables a string vector to specify the column names of independent variables included in the regression model other than the variableOfInterst.
#' @param studyID a string to specify the column name of study ID.
#' @param data a data frame containing the data
#' @param CILevel a number between 0 to 1, used to specify the confidence interval to be plotted in the forrest plot.
#' @param main a string to specify the title of the forrest plot
#' @param ifScale a vector of two logic values, specifing if the dependent variable and the variableOfInterst should be scaled when calculating the effect size.
#' @param cex a number specifying the amount by which plotting text and symbols should be scaled relative to the default in the forrest plot.
#' @return returns data frame describing the effect size of variableOfInterst on value in each individule studies, as well as the over all effect size.
#' @export
metaAnalysis=function(value,variableOfInterst,otherVariables,
                      studyID,data,CILevel,main,
                      ifScale=c(T,F),cex=1){
  study_result=NULL
  for(std in unique(data[,studyID])){
    sub_data=subset(data,data[,studyID]==std)
    sub_data=na.omit(sub_data)
    NL=apply(sub_data[,c(value,variableOfInterst,otherVariables)],2,function(x){length(unique(x))})
    if(!all(NL>1)){cat(std,"is skipped. One of the variable have 0 variance.\n");next}
    if(ifScale[1]){sub_data[,value]=scale(sub_data[,value])}
    if(ifScale[2]){sub_data[,variableOfInterst]=scale(sub_data[,variableOfInterst])}
    x=sub_data[,c(value,variableOfInterst,otherVariables)]
    colnames(x)[1]="Y"
    LM=lm(Y~.,data=x)
    CE=t(summary(LM)$coefficients[2,])
    CI=confint(LM,parm=2,level=CILevel)
    t1=cbind("study_id"=std,data.frame(CE,check.names=F),
             data.frame(CI,check.names=F),"N"=nrow(sub_data))
    study_result=rbind(study_result,t1)
  }
  study_result=na.omit(study_result)
  res=metafor::rma.uni(yi=study_result$Estimate, vi=(study_result$`Std. Error`)^2)
  metafor::forest(res, slab=study_result$study_id,main=main,
                  xlab="Effect Size", mlab="RE Model for All Studies",cex=cex)
  t1=data.frame("Summary",res$b[1],res$se,res$zval,res$pval,res$ci.lb,res$ci.ub,
                sum(study_result$N))
  names(t1)=names(study_result)
  study_result=rbind(study_result,t1)
  rownames(study_result)=NULL
  return(study_result)
}

#' perform generalized linear model analysis to estimate effect size.
#'
#' A function that perform generalized linear model analysis to estimate effect size.
#' @param value a string to specify the column name of the dependent variable (y)
#' @param variableOfInterst a string to specify the column name of the independent variable of interest (x1)
#' @param otherVariables a string vector to specify the column names of independent variables included in the regression model other than the variableOfInterst.
#' @param parameter a string to specify what summary statistics is the dependent variale.
#' @param studyID a string to specify the column name of study ID.
#' @param label a string to specify the name the column that contains the cluster label or name.
#' @param data a data frame containing the data. Usually a long form data frame returned by collectData.
#' @param CILevel a number between 0 to 1, used to specify the confidence interval to be plotted in the forrest plot.
#' @param ifScale a vector of two logic values, specifing if the dependent variable and the variableOfInterst should be scaled when calculating the effect size.
#' @return returns data frame describing the overall effect size of variableOfInterst on value.  May be slightly different from the value reported from the function metaAnalysis.
#' @details The function use the model value ~ variableOfInterst + otherVariables + studyID to estimate the effect size. Use it as a screening tool. Use metaAnalysis function to analyze an  effect size in more detail.
#' @export
glmAnalysis=function(value="value",variableOfInterst="SUBJECT_AGE",parameter,
                     otherVariables=c("GENDER"),studyID="study",label="label",
                     data,CILevel=0.95,ifScale=c(T,F)){
  result=NULL
  for(L in unique(data[,label])){
    sub_data=data[data$parameter_name==parameter&data[,label]==L,]
    sub_data=na.omit(sub_data)
    if(length(unique(sub_data[,studyID]))>1){
      sub_data=sub_data[,c(value,variableOfInterst,otherVariables,studyID)]
    }else{sub_data=sub_data[,c(value,variableOfInterst,otherVariables)]}

    NL=apply(sub_data,2,function(x){length(unique(x))})
    if(!all(NL>1)){cat(L,"is skipped. One of the variable have 0 variance.\n");next}
    if(ifScale[1]){sub_data[,value]=scale(sub_data[,value])}
    if(ifScale[2]){sub_data[,variableOfInterst]=scale(sub_data[,variableOfInterst])}
    colnames(sub_data)[1]="Y"
    LM=lm(Y~.,data=sub_data)
    CE=t(summary(LM)$coefficients[2,])
    CI=confint(LM,parm=2,level=CILevel)
    t1=cbind("label"=L,data.frame(CE,check.names=F),
             data.frame(CI,check.names=F),"N"=nrow(sub_data))
    result=rbind(result,t1)
  }
  rownames(result)=NULL
  colnames(result)=c("label","Effect_size","SE","t_value","p_value","lower","upper","N")
  return(result)
}


#' plot the result from the glmAnalysis function
#'
#' A function that plot the result from the glmAnalysis function.
#' @param GA a data frame returned from the function glmAnalysis.
#' @param size font size of texts in the plot
#' @return NULL
#' @export
plotGA=function(GA,size=16){
  GA$label=factor(GA$label,levels=GA$label)
  p <- ggplot2::ggplot(GA, ggplot2::aes(y=GA$Effect_size,x=GA$label,ymin=lower, ymax=upper))+
    ggplot2::geom_pointrange()+
    ggplot2::geom_hline(yintercept = 0, linetype=2)+
    ggplot2::coord_flip()+
    ggplot2::xlab('label')+
    ggplot2::ylab('Effect size')+
    ggplot2::theme(text=ggplot2::element_text(size=size))
  print(p)
}


# Define clustering functions---------

#' cluster cytometry data using hierarchical clustering
#'
#' A function that cluster cytometry data using hierarchical clustering.
#' @param fcsFrame a flow frame.
#' @param excludeClusterParameters A vector specifiying the name of markers not to be used for clustering.
#' @param minimumClusterSizePercent a number between 0 and 1, used to specify the minimum size of a cluster relative to all events.
#' @return a list of clusters. Each cluster contains the ID of all cells that belong to the cluster.
#' @export
flowHC=function(fcsFrame,excludeClusterParameters,minimumClusterSizePercent=0.05){
  antibodies=markerFinder(fcsFrame)
  excludeClusterParameters=toupper(excludeClusterParameters)
  excludeClusterParameters=union(excludeClusterParameters,c("TIME","SAMPLE_ID"))
  w=which(!antibodies%in%excludeClusterParameters)
  expr=flowCore::exprs(fcsFrame)
  expr=expr[,w,drop=F]
  D=dist(expr)
  HC=fastcluster::hclust(D,method="ward.D")
  CL=HC.assign(HC,minimumClusterSizePercent)
  return(CL)
}


#' cluster cytometry data using FlowSOM
#'
#' A function that cluster cytometry data using FlowSOM.
#' @param fcsFrame a flow frame.
#' @param excludeClusterParameters A vector specifiying the name of markers not to be used for clustering.
#' @return a list of clusters. Each cluster contains the ID of all cells that belong to the cluster.
#' @export
flowSOM.MC=function(fcsFrame,excludeClusterParameters){
  antibodies=markerFinder(fcsFrame)
  excludeClusterParameters=toupper(excludeClusterParameters)
  excludeClusterParameters=union(excludeClusterParameters,c("TIME","SAMPLE_ID"))
  w=which(!antibodies%in%excludeClusterParameters)
  fSOM <- FlowSOM::ReadInput(fcsFrame, transform = FALSE, scale = FALSE)
  fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse = w,
                            xdim = 10, ydim = 10)
  fSOM <- FlowSOM::BuildMST(fSOM)
  meta <- FlowSOM::metaClustering_consensus(fSOM$map$codes,k=40)
  CL <- meta[fSOM$map$mapping[, 1]]
  CL=lapply(unique(CL),function(x){which(CL==x)})
  return(CL)
}


# Define back-end functions----------
findCutoff=function(x,returnSil=F,useBL=T,minX=0){
  if(length(x)>1000){x=sample(x,1000)}
  if(useBL==T){x=baselineCut(x)}
  #valley=sapply(seq(0.01,0.99,length.out=100),function(q){quantile(x,q)})
  valley=seq(min(x),max(x),length.out=100)
  if(!is.null(minX)){valley=valley[valley>minX]}
  D=dist(x)
  sil=sapply(valley,function(v){
    cluster=1*(x>v)+1
    if(length(unique(cluster))<2){return(-2)}
    ss <- cluster::silhouette(cluster,D)
    return(mean(ss[, 3]))
  })
  valley=valley[which.max(sil)]
  if(returnSil==F){return(valley)}else{return(c("cuoff"=valley,"sil"=max(sil)))}
}
labelUnifier=function(clusterLabel){
  clusterLabel=strsplit(as.character(clusterLabel),split="\\|")
  clusterLabel=lapply(clusterLabel,sort)
  clusterLabel=sapply(clusterLabel,paste,collapse="|")
  clusterLabel=toupper(clusterLabel)
  return(clusterLabel)
}
labelCombiner=function(labels){
  label1=labels[1]
  label2=labels[2]
  if(nchar(label1)!=nchar(label2)){return(NA)}
  m1=gsub("\\+|-","",label1)
  m2=gsub("\\+|-","",label2)
  if(m1!=m2){return(NA)}
  label1=strsplit(label1,"\\|")[[1]]
  label2=strsplit(label2,"\\|")[[1]]
  s1=sapply(label1, function(x){substring(x,nchar(x))})
  s2=sapply(label2, function(x){substring(x,nchar(x))})

  if(sum(s1!=s2)!=1){return(NA)}
  w=which(s1==s2)
  label=label1[w]
  label=paste(label,collapse="|")
  return(label)
}
HC.assign=function(HC,minimumClusterSizePercent){
  CL=list()
  j=0
  NL=min(100,nrow(HC$merge))
  for(k in 2:100){
    t1=cutree(HC,k)
    w=lapply(unique(t1),function(x){which(t1==x)})
    w=w[sapply(w, length)>minimumClusterSizePercent*(nrow(HC$merge)+1)]

    if(length(w)==0){break}
    CL_new=c(CL,w)
    if(length(unique(CL_new))==length(unique(CL))){j=j+1}else{j=0}
    CL=unique(CL_new)
    if(j>20){break}
  }
  return(CL)
}
baselineCut=function(x){
  breaks=seq(min(x),max(x),length.out=100)
  h=hist(x,breaks,plot=F)
  #h2=h$counts[h$counts>1]
  baseline=max(3,mean(h$counts)/10)
  bound=range(which(h$counts>baseline))
  x=x[x>breaks[bound[1]]&x<breaks[bound[2]]]
  return(x)
}
trisect=function(x){
  if(length(x)>2000){x=sample(x,2000)}
  x=baselineCut(x)
  D=dist(x)
  cutoff1=findCutoff(x)

  t1=findCutoff(x[x>cutoff1],useBL=F,minX=NULL)
  cluster=.bincode(x,sort(c(-Inf,cutoff1,t1,+Inf),decreasing=F))
  ss <- cluster::silhouette(cluster,D)
  s1=mean(ss[, 3])

  t2=findCutoff(x[x<cutoff1],useBL=F,minX=NULL)
  cluster=.bincode(x,sort(c(-Inf,cutoff1,t2,+Inf),decreasing=F))
  ss <- cluster::silhouette(cluster,D)
  s2=mean(ss[, 3])

  if(s2>s1){cutoff2=t2}else{cutoff2=t1}
  return(sort(c(cutoff1,cutoff2),decreasing=T))
}


# Archive functions------
labelSummary=function(allData,minStudy=2){
  result=NULL
  for(label in unique(allData$label)){
    mid_data=allData[allData$label==label,]
    for(param in unique(allData$parameter_name)){
      sub_data=mid_data[mid_data$parameter_name==param,]
      if(nrow(sub_data)<1){next}
      N_study=length(unique(sub_data$study))
      if(N_study<minStudy){next}
      t1=data.frame("label"=label,"parameter_name"= param,"number_of_study"=N_study,"sample_size"=nrow(sub_data),
                    "studies"=paste(unique(sub_data$study),collapse="|"), row.names=NULL )
      result=rbind(result,t1)
    }
  }
  return(result)
}
