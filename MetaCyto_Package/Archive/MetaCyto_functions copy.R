# Define user level functions----------------
#' pre-processing fcs files by transformation and compensation.
#'
#' It transform and compensate for the raw fcs files and write out the processed data to a new set of fcs files.
#' @param inputMeta A dataframe containing 2 columns: a column called "fcs_files" that contains the location of each fcs file on the hard drive and a column called "study_id" that specify what study each fcs file belongs to.
#' @param assay Either "FCM" or "CyTOF" to indicate the type of cytometry data.
#' @param outpath a string indicating the directory the pre-processed fcs files will be wrote to.
#' @param b a positive number used to specify the arcsinh transformtion. f(x) = asinh (b*x) where x is the original value and f(x) is the value after transformation. The suggested value is 1/150 for flow cytometry (FCM) data and 1/8 for CyTOF data.
#' @param fileSampleSize An integer specifying the number of events sampled from each fcs file. If NULL, all the events will be pre-processed and wrote out to the new fcs files.
#' @param excludeTransformParameters A vector specifiying the name of parameters not to be transformed (left at linear scale).
#' @return Does not return anything. The output is wrote to the directory specified by the "outpath".
#' @details The function takes a data frame which specify the location of the fcs files and the panels the fcs files belong to. It transform the cytometry data using the arcsinh transformtion. For flow cytometry data, it compensate the data using the compensation matrix supplied in the fcs file. the preprocessed fcs files and a table called "processed_sample_summary.csv" will be wrote out to outpath as well.
#' @export

preprocessing = function(inputMeta,
                         assay=c("FCM", "CyTOF"),
                         outpath,
                         b=1/200,
                         fileSampleSize=5000,
                         excludeTransformParameters=c("FSC-A","FSC-W","FSC-H","Time","Cell_length")){

  fcs_param=NULL
  fcs_names= gsub(".*/","",inputMeta$fcs_files)
  #fcs_name=inputMeta$fcs_files
  inputMeta=cbind(fcs_names,inputMeta)
  excludeTransformParameters=paste(excludeTransformParameters,collapse="|")

  for(std in unique(inputMeta$study_id)){

    cat("Preprocessing, study ID = ",std, "\n")

    # 1) identify sample :
    fcs_files=subset(inputMeta$fcs_files,inputMeta$study_id==std)
    fcs_files=as.character(fcs_files)

    # 2) read the fcs files
    fcs = flowCore::read.flowSet(fcs_files,transformation="linearize",alter.names=FALSE)
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

    # 5) outputting resuts
    flowCore::write.flowSet(fcs, outdir = paste0(outpath,'/',std))

    #a function to extract parameter names
    fp = function(frame) {
      # get channel names
      channels=as.vector(flowCore::pData(flowCore::parameters(frame))[,1])
      # get antibody names
      antibodies=as.vector(flowCore::pData(flowCore::parameters(frame))[,2])
      antibodies=sapply(1:length(antibodies), function(i){
        if(is.na(antibodies[i])|antibodies[i]=="NA"){return(channels[i])}else{antibodies[i]}
      })

      channels = paste(channels,collapse='|' )
      antibodies = paste(antibodies,collapse='|' )
      return(c(channels,antibodies))
    }
    param= flowCore::fsApply(fcs, fp)
    param= toupper(param)
    fcs_param=rbind(fcs_param,data.frame("fcs_names"=rownames(param),"channels"=param[,1],"antibodies"=param[,2]))
  }#end of each std
  inputMeta=dplyr::left_join(fcs_param,inputMeta,by="fcs_names")
  #inputMeta=select(inputMeta,-fcs_files)
  write.csv(inputMeta,paste0(outpath,"/processed_sample_summary.csv"),row.names=F)
  cat("Preprocess result stored in the folder:",outpath, "\n")
}#end of function

#' @export
autoCluster= function(preprocessOutputFolder,
                      outpath="cluster_output",
                      excludeClusterParameters=c("TIME"),
                      minimumClusterSizePercent=0.05,
                      fileSampleSize=1000,
                      labelQuantile=0.8){
  #read the output from preprocessing
  inputMeta=read.csv(paste0(preprocessOutputFolder,'/processed_sample_summary.csv'),stringsAsFactors=F)
  #create output foler
  dir.create(file.path(outpath),recursive=T)

  #prepare exclude parameters
  excludeClusterParameters=toupper(excludeClusterParameters)
  cluster_summary=NULL
  for(std in unique(inputMeta$study_id)){
    cat("Clustering , study ID = ",std, "\n")

    ##### 1) read sample files for each study--------
    fcs_files=subset(inputMeta$fcs_names,inputMeta$study_id==std)
    fcs_files=paste(preprocessOutputFolder,std,fcs_files,sep="/")
    fcs=flowCore::read.flowSet(fcs_files, transformation = "linearize", alter.names = FALSE)
    # make sure the fcs file antibody names are the same as the preprocessed output
    antibodies=subset(inputMeta$antibodies,inputMeta$study_id==std)[1]
    antibodies=strsplit(antibodies,"\\|")[[1]]

    ##### 2) subset the cells in fcs-----------
    # Get expression matrix
    expr=flowCore::fsApply(fcs,function(x){
      v=flowCore::exprs(x);
      if(nrow(v)>fileSampleSize){
        v=v[sample(1:nrow(v),fileSampleSize,replace=T),]
      }
      return(v)
    })
    colnames(expr)=antibodies
    expr_sample=ceiling(seq_along(1:nrow(expr))/fileSampleSize)


    # subset on columns
    w=!antibodies%in%excludeClusterParameters
    antibodies=antibodies[w]
    if(length(antibodies)<2){next}
    expr=expr[,w]
    expr_scale=scale(expr,center=F,scale=T)

    # 3) cluster the subset samples, find clusters large enough------------
    max_cell_N=20000
    if(nrow(expr_scale)>max_cell_N){
      ind=sample(1:nrow(expr_scale),max_cell_N,replace=F)

      D=dist(expr_scale[ind,])
    }else{
      D=dist(expr_scale)
    }
    HC=fastcluster::hclust(D,method="ward.D")
    CL=HC.assign(HC,minimumClusterSizePercent)

    # 4) label the clusters:--------------
    #find cutoff of each parameter
    cutoff=apply(X=expr,MARGIN=2,FUN=findCutoff)

    #label each cluster
    CL_label=rep(NA,length(CL))
    if(nrow(expr_scale)>max_cell_N){expr_ind=expr[ind,]}else{expr_ind=expr}
    no_label=sapply(1:length(antibodies),function(i){
      x1=quantile(expr_ind[,i],minimumClusterSizePercent)
      x2= quantile(expr_ind[,i],(1-minimumClusterSizePercent))
      if(cutoff[i]<x1|cutoff[i]>x2){return(T)}else{return(F)}
    })
    for(i in 1:length(CL)){
      l=c()
      for(j in 1:length(antibodies)){
        if(no_label[j]==T){next}
        x=expr_ind[CL[[i]],j]
        if(quantile(x,labelQuantile)<cutoff[j]){
          l=c(l,paste0(antibodies[j],"-"))
        }else if(quantile(x,(1-labelQuantile))>cutoff[j]){
          l=c(l,paste0(antibodies[j],"+"))
        }
      }
      l=paste(l,collapse="|")
      CL_label[i]=l
    }

    ##### 5) merge clusters-----------
    #merge clusters with the same labels
    temp=list()
    for(l in unique(CL_label)){
      w=which(CL_label==l)
      temp= c(temp, list(unique(unlist(CL[w]))) )
    }
    CL=temp
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
      new_CL=lapply(w, function(i){
        l1=label_pair[1,i]
        l2=label_pair[2,i]
        w1=which(CL_label==l1)
        w2=which(CL_label==l2)
        new_CL=unique(unlist(CL[c(w1,w2)]))
        return(new_CL)
      })
      CL=c(CL,new_CL)
    }
    #make sure CL is unique
    temp=list()
    for(l in unique(CL_label)){
      w=which(CL_label==l)
      temp= c(temp, list(unique(unlist(CL[w]))) )
    }
    CL_label=unique(CL_label)
    CL_label=labelUnifier(CL_label)
    CL=temp

    # 6) assign cluster to all cells-------------
    CL_label=CL_label[CL_label!=""]
    if(length(CL_label)==0){print("No clusters found");next}
    CL_label2=strsplit(CL_label,split="\\|")
    CL=list()
    for(i in 1:length(CL_label2)){
      rows=sapply(CL_label2[[i]],function(x){
        m=gsub("\\+$|-$","",x)
        if(grepl("\\+$",x)){w=which(expr[,m]>cutoff[m])}
        if(grepl("-$",x)){w=which(expr[,m]<cutoff[m])}
        #w=which(expr[,m]<quantile(expr_ind[CL[[i]],m],1) & expr[,m]>quantile(expr_ind[CL[[i]],m],0) )
        return(w)
      })
      if(length(CL_label2[[i]])>1){rows=Reduce(intersect,rows)}
      CL[[i]]=rows
    }
    id_keep=which(sapply(CL,length)>minimumClusterSizePercent*nrow(expr))
    if(length(id_keep)<1){next}
    CL=CL[id_keep]
    CL_label=CL_label[id_keep]

    ### 6) record cluster data-------
    dir.create(file.path(paste(outpath,std,sep="/")))
    #find cluster median
    CL_median=lapply(CL,function(x){apply(expr[x,],2,median)})
    CL_median=data.frame(matrix(unlist(CL_median), ncol=length(CL_median[[1]]), byrow=T))
    CL_median=cbind(std,paste0("cluster",1:length(CL)),unique(CL_label),CL_median)
    names(CL_median)=c("study","cluster_id","label",antibodies)
    write.csv(CL_median,paste(outpath,std,"cluster_median.csv",sep="/"),row.names=F)
    cluster_summary=rbind(cluster_summary,data.frame("Study"=std,"cluster_id"=CL_median$cluster_id,"label"=unique(CL_label)))

    ### 7) record sample data-----------
    sample_name=subset(inputMeta$fcs_files,inputMeta$study_id==std)
    SP_stat=NULL
    for(i in 1:length(CL)){
      for(j in unique(expr_sample)){
        w=intersect( which(expr_sample==j), CL[[i]] )
        frac=length(w)/fileSampleSize
        if(length(w)>1){M=apply(expr[w,],2,median)}
        if(length(w)==1){M=expr[w,]}
        if(length(w)==0){M=rep(NA,ncol(expr))}
        t1=c(std,paste0("cluster",i),CL_label[i],sample_name[j],M,frac)
        SP_stat=rbind(SP_stat,t1)
      }
    }
    colnames(SP_stat)=c("study","cluster_id","label","fcs_files",antibodies,"fraction")
    write.csv(SP_stat,paste(outpath,std,"cluster_stats_in_each_sample.csv",sep="/"),row.names=F)


    ## 8) plot the density plot--------
    filename=paste(outpath,std,"density_plot.pdf",sep="/")
    height=3*length(CL);width=3*length(antibodies)
    pdf(filename,width=width,height=height)
    par( mfcol = c(length(CL), length(antibodies) ) )
    for(i in 1:length(antibodies)){
      x_all=expr[,i]
      for(j in 1:length(CL)){
        b=seq(min(x_all),max(x_all), ((max(x_all)-min(x_all))/100) )
        hist(x_all,col=rgb(0, 0, 0, 0.2),xlab=antibodies[i],breaks=b,freq=T,border=F,main=paste0(antibodies[i]," of cluster ",j," : \n", CL_label[j]))
        hist(expr[CL[[j]],i],add=T,breaks=b,col=rgb(1, 0, 0, 0.5),freq=T,border=F)
        abline(v=cutoff[i])
      }
    }
    dev.off()
  }#end of each study
  #write.csv(cluster_summary,paste(outpath,"cluster_list.csv",sep="/"),row.names=F)
}

#' @export
searchCluster=function(preprocessOutputFolder,
                       outpath="search_output",
                       clusterLabel,
                       fileSampleSize=1000,
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

    ##### 1) read sample files for each study############################################################
    fcs_files=subset(inputMeta$fcs_names,inputMeta$study_id==std)
    fcs_files=paste(preprocessOutputFolder,std,fcs_files,sep="/")
    fcs=flowCore::read.flowSet(fcs_files, transformation = "linearize", alter.names = FALSE)
    # make sure the fcs file antibody names are the same as the preprocessed output
    antibodies=subset(inputMeta$antibodies,inputMeta$study_id==std)[1]
    antibodies=strsplit(antibodies,"\\|")[[1]]

    ##### 2) subset the cells in fcs ############################################################
    # Get expression matrix
    expr=flowCore::fsApply(fcs,function(x){
      v=flowCore::exprs(x);
      if(nrow(v)>fileSampleSize){
        v=v[sample(1:nrow(v),fileSampleSize,replace=T),]
      }
      return(v)
    })
    colnames(expr)=antibodies
    expr_sample=ceiling(seq_along(1:nrow(expr))/fileSampleSize)

    ##### 3) bisect each marker############################################################
    #find cluster whose markers are included in the study
    CL_label=clusterLabel
    CL_label2=strsplit(clusterLabel,split="&|\\|")
    in_study=sapply(CL_label2,function(x){
      Ab_in_label=gsub("\\+$|-$|\\^NE$|\\^LO$|\\^HI$","",x)
      if(length(setdiff(Ab_in_label,antibodies))>0){return(F)}else{return(T)}
    })
    CL_label=CL_label[in_study]
    CL_label2=CL_label2[in_study]
    if(length(CL_label)<1){next}
    #find if trisect is needed
    t1=unlist(CL_label2)
    t1=t1[grep("\\^",t1)]
    t1=gsub("\\^NE$|\\^LO$|\\^HI$","",t1)
    Ab_tri=unique(t1)

    #find cutoff of each parameter
    cutoff=apply(X=expr,MARGIN=2,FUN=findCutoff)
    if(length(Ab_tri)>0){triS=apply(expr,2,trisect)}

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
    CL_count=sapply(CL,length)
    CL=CL[CL_count>1]
    CL_label=CL_label[CL_count>1]
    if(length(CL)<1){next}

    ##### 4) record cluster data#########################################
    dir.create(file.path(paste(outpath,std,sep="/")))
    #find cluster median
    CL_median=lapply(CL,function(x){apply(expr[x,],2,median)})
    CL_median=data.frame(matrix(unlist(CL_median), ncol=length(CL_median[[1]]), byrow=T))
    CL_median=cbind(std,paste0("cluster",1:length(CL)),unique(CL_label),CL_median)
    names(CL_median)=c("study","cluster_id","label",antibodies)
    write.csv(CL_median,paste(outpath,std,"cluster_median.csv",sep="/"),row.names=F)
    cluster_summary=rbind(cluster_summary,data.frame("Study"=std,"cluster_id"=CL_median$cluster_id,"label"=unique(CL_label)))

    ##### 5) record sample data#########################################
    sample_name=subset(inputMeta$fcs_files,inputMeta$study_id==std)
    SP_stat=NULL
    for(i in 1:length(CL)){
      for(j in unique(expr_sample)){
        w=intersect( which(expr_sample==j), CL[[i]] )
        frac=length(w)/fileSampleSize
        if(length(w)>1){M=apply(expr[w,],2,median)}
        if(length(w)==1){M=expr[w,]}
        if(length(w)==0){M=rep(NA,ncol(expr))}
        t1=c(std,paste0("cluster",i),CL_label[i],sample_name[j],M,frac)
        SP_stat=rbind(SP_stat,t1)
      }
    }
    colnames(SP_stat)=c("study","cluster_id","label","fcs_files",antibodies,"fraction")
    write.csv(SP_stat,paste(outpath,std,"cluster_stats_in_each_sample.csv",sep="/"),row.names=F)


    ##### 6) plot the density plot############################
    if(ifPlot==T){
      filename=paste(outpath,std,"density_plot.pdf",sep="/")
      height=3*length(CL);width=3*length(antibodies)
      pdf(filename,width=width,height=height)
      par( mfcol = c(length(CL), length(antibodies) ) )
      for(i in 1:length(antibodies)){
        x_all=expr[,i]
        for(j in 1:length(CL)){
          b=seq(min(x_all),max(x_all), ((max(x_all)-min(x_all))/100) )
          hist(x_all,col=rgb(0, 0, 0, 0.2),xlab=antibodies[i],breaks=b,freq=T,border=F,main=paste0(antibodies[i]," of cluster ",j," : \n", CL_label[j]))
          hist(expr[CL[[j]],i],add=T,breaks=b,col=rgb(1, 0, 0, 0.5),freq=T,border=F)
          abline(v=cutoff[i])
          if(antibodies[i]%in%Ab_tri){abline(v=triS[,i])}
        }
      }
      dev.off()
    }
  }#end of each study
  #write.csv(cluster_summary,paste(outpath,"cluster_list.csv",sep="/"),row.names=F)
}
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
      fcs=flowCore::read.FCS(x)
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

#' @export
collectData=function(files,longform=T){
  all_data=NULL
  nms=NULL
  for(fn in files){
    #read output data from "autoCluster" function
    cluster_stat=read.csv(fn,stringsAsFactors=F,check.names=F)
    if(is.null(nms)){nms=colnames(cluster_stat)}
    w=which(nms=="fcs_files")+1
    if(longform==T){
      cluster_stat=tidyr::gather(cluster_stat,parameter_name, value ,w:ncol(cluster_stat))
    }else(colnames(cluster_stat)=nms)
    all_data=rbind(all_data,cluster_stat)
  }
  return(all_data)
}

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

#' @export
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
# Define back-end functions----------
#' @export
findCutoff=function(x,returnSil=F,useBL=T,minX=0){
  if(length(x)>2000){x=sample(x,2000)}
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




# Experimental functions----------

#'
#' @export
autoCluster.oneSample= function(fcsFile,
                      outpath="cluster_output",
                      excludeClusterParameters=c("TIME"),
                      minimumClusterSizePercent=0.05,
                      labelQuantile=0.8,
                      clusterFunction=flowHC,...){

  #create output foler
  dir.create(file.path(outpath),recursive=T)

  #prepare exclude parameters
  excludeClusterParameters=toupper(excludeClusterParameters)
  cluster_summary=NULL

    cat("Clustering , fcs file name = ",fcsFile, "\n")

    ##### 1) read sample files for each study--------
    fcs=flowCore::read.FCS(fcsFile, transformation = "linearize", alter.names = FALSE,truncate_max_range=F)
    # make sure the fcs file antibody names are the same as the preprocessed output
    channels=as.vector(flowCore::pData(flowCore::parameters(fcs))[,1])
    # get antibody names
    antibodies=as.vector(flowCore::pData(flowCore::parameters(fcs))[,2])
    antibodies=sapply(1:length(antibodies), function(i){
      if(is.na(antibodies[i])|antibodies[i]=="NA"){return(channels[i])}else{antibodies[i]}
    })
    antibodies=toupper(antibodies)

    ##### 2) subset the cells in fcs-----------
    # Get expression matrix
    expr=flowCore::exprs(fcs);
    colnames(expr)=antibodies

    # subset on columns
    w=!antibodies%in%excludeClusterParameters
    antibodies=antibodies[w]
    expr=expr[,w]
    expr_scale=scale(expr,center=F,scale=T)
    expr_sample=rep(1,nrow(expr))

    # 3) cluster the subset samples, find clusters large enough------------
    
    max_cell_N=20000
    if(nrow(expr_scale)>max_cell_N){
       ind=sample(1:nrow(expr_scale),max_cell_N,replace=F)
     }else{
       ind=1:nrow(expr_scale)
     }
    CL=clusterFunction(exprM=expr_scale[ind,],...)
    
    # HC=fastcluster::hclust(D,method="ward.D")
    # CL=HC.assign(HC,minimumClusterSizePercent)

    # 4) label the clusters:--------------
    #find cutoff of each parameter
    cutoff=apply(X=expr,MARGIN=2,FUN=findCutoff)

    #label each cluster
    CL_label=rep(NA,length(CL))
    if(nrow(expr_scale)>max_cell_N){expr_ind=expr[ind,]}else{expr_ind=expr}
    no_label=sapply(1:length(antibodies),function(i){
      x1=quantile(expr_ind[,i],minimumClusterSizePercent)
      x2= quantile(expr_ind[,i],(1-minimumClusterSizePercent))
      if(cutoff[i]<x1|cutoff[i]>x2){return(T)}else{return(F)}
    })
    for(i in 1:length(CL)){
      l=c()
      for(j in 1:length(antibodies)){
        if(no_label[j]==T){next}
        x=expr_ind[CL[[i]],j]
        if(quantile(x,labelQuantile)<cutoff[j]){
          l=c(l,paste0(antibodies[j],"-"))
        }else if(quantile(x,(1-labelQuantile))>cutoff[j]){
          l=c(l,paste0(antibodies[j],"+"))
        }
      }
      l=paste(l,collapse="|")
      CL_label[i]=l
    }

    ##### 5) merge clusters-----------
    #merge clusters with the same labels
    temp=list()
    for(l in unique(CL_label)){
      w=which(CL_label==l)
      temp= c(temp, list(unique(unlist(CL[w]))) )
    }
    CL=temp
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
      new_CL=lapply(w, function(i){
        l1=label_pair[1,i]
        l2=label_pair[2,i]
        w1=which(CL_label==l1)
        w2=which(CL_label==l2)
        new_CL=unique(unlist(CL[c(w1,w2)]))
        return(new_CL)
      })
      CL=c(CL,new_CL)
    }
    #make sure CL is unique
    temp=list()
    for(l in unique(CL_label)){
      w=which(CL_label==l)
      temp= c(temp, list(unique(unlist(CL[w]))) )
    }
    CL_label=unique(CL_label)
    CL_label=labelUnifier(CL_label)
    CL=temp

    # 6) assign cluster to all cells-------------
    CL_label=CL_label[CL_label!=""]
    CL_label2=strsplit(CL_label,split="\\|")
    CL=list()
    for(i in 1:length(CL_label2)){
      rows=sapply(CL_label2[[i]],function(x){
        m=gsub("\\+$|-$","",x)
        if(grepl("\\+$",x)){w=which(expr[,m]>cutoff[m])}
        if(grepl("-$",x)){w=which(expr[,m]<cutoff[m])}
        #w=which(expr[,m]<quantile(expr_ind[CL[[i]],m],1) & expr[,m]>quantile(expr_ind[CL[[i]],m],0) )
        return(w)
      })
      if(length(CL_label2[[i]])>1){rows=Reduce(intersect,rows)}
      CL[[i]]=rows
    }
    id_keep=which(sapply(CL,length)>minimumClusterSizePercent*nrow(expr))
    if(length(id_keep)<1){next}
    CL=CL[id_keep]
    CL_label=CL_label[id_keep]

    ### 6) record cluster data-------
    #find cluster median
    CL_median=lapply(CL,function(x){apply(expr[x,],2,median)})
    CL_median=data.frame(matrix(unlist(CL_median), ncol=length(CL_median[[1]]), byrow=T))
    CL_median=cbind(paste0("cluster",1:length(CL)),unique(CL_label),CL_median)
    names(CL_median)=c("cluster_id","label",antibodies)
    write.csv(CL_median,paste(outpath,"cluster_median.csv",sep="/"),row.names=F)

    ### 7) record sample data-----------
    SP_stat=NULL
    for(i in 1:length(CL)){
      for(j in unique(expr_sample)){
        w=intersect( which(expr_sample==j), CL[[i]] )
        frac=length(w)/nrow(expr)
        if(length(w)>1){M=apply(expr[w,],2,median)}
        if(length(w)==1){M=expr[w,]}
        if(length(w)==0){M=rep(NA,ncol(expr))}
        t1=c(paste0("cluster",i),CL_label[i],fcsFile,M,frac)
        SP_stat=rbind(SP_stat,t1)
      }
    }
    colnames(SP_stat)=c("cluster_id","label","fcs_files",antibodies,"fraction")
    write.csv(SP_stat,paste(outpath,"cluster_stats_in_each_sample.csv",sep="/"),row.names=F)


    ## 8) plot the density plot--------
    filename=paste(outpath,"density_plot.pdf",sep="/")
    height=3*length(CL);width=3*length(antibodies)
    pdf(filename,width=width,height=height)
    par( mfcol = c(length(CL), length(antibodies) ) )
    for(i in 1:length(antibodies)){
      x_all=expr[,i]
      for(j in 1:length(CL)){
        b=seq(min(x_all),max(x_all), ((max(x_all)-min(x_all))/100) )
        hist(x_all,col=rgb(0, 0, 0, 0.2),xlab=antibodies[i],breaks=b,freq=T,border=F,main=paste0(antibodies[i]," of cluster ",j," : \n", CL_label[j]))
        hist(expr[CL[[j]],i],add=T,breaks=b,col=rgb(1, 0, 0, 0.5),freq=T,border=F)
        abline(v=cutoff[i])
      }
    }
    dev.off()

    ## 9) write.out cell id that belong to each clusters--------
    filename=paste(outpath,"cluster_cell_id.RData",sep="/")
    names(CL)=CL_label
    save(CL,file=filename)

  #write.csv(cluster_summary,paste(outpath,"cluster_list.csv",sep="/"),row.names=F)
}
#'
#' @export
searchCluster.oneSample=function(fcsFile,
                       outpath="search_output",
                       clusterLabel,
                       ifPlot=T){

  #create output foler
  dir.create(file.path(outpath),recursive=T)

  #prepare exclude parameters
  clusterLabel=toupper(clusterLabel)

    cat("Searching , file name = ",fcsFile, "\n")

    ##### 1) read sample files for each study############################################################
    fcs=flowCore::read.FCS(fcsFile, transformation = "linearize", alter.names = FALSE,truncate_max_range=F)
    # make sure the fcs file antibody names are the same as the preprocessed output
    channels=as.vector(flowCore::pData(flowCore::parameters(fcs))[,1])
    # get antibody names
    antibodies=as.vector(flowCore::pData(flowCore::parameters(fcs))[,2])
    antibodies=sapply(1:length(antibodies), function(i){
      if(is.na(antibodies[i])|antibodies[i]=="NA"){return(channels[i])}else{antibodies[i]}
    })
    antibodies=toupper(antibodies)

    ##### 2) subset the cells in fcs ############################################################
    # Get expression matrix
    expr=flowCore::exprs(fcs);
    colnames(expr)=antibodies
    expr_sample=rep(1,nrow(expr))

    ##### 3) bisect each marker############################################################
    #find cluster whose markers are included in the study
    CL_label=clusterLabel
    CL_label2=strsplit(clusterLabel,split="&|\\|")
    in_study=sapply(CL_label2,function(x){
      Ab_in_label=gsub("\\+$|-$|\\^NE$|\\^LO$|\\^HI$","",x)
      if(length(setdiff(Ab_in_label,antibodies))>0){return(F)}else{return(T)}
    })
    CL_label=CL_label[in_study]
    CL_label2=CL_label2[in_study]
    if(length(CL_label)<1){next}
    #find if trisect is needed
    t1=unlist(CL_label2)
    t1=t1[grep("\\^",t1)]
    t1=gsub("\\^NE$|\\^LO$|\\^HI$","",t1)
    Ab_tri=unique(t1)

    #find cutoff of each parameter
    cutoff=apply(X=expr,MARGIN=2,FUN=findCutoff)
    if(length(Ab_tri)>0){triS=apply(expr,2,trisect)}

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
    CL_count=sapply(CL,length)
    CL=CL[CL_count>1]
    CL_label=CL_label[CL_count>1]
    if(length(CL)<1){next}

    ##### 4) record cluster data#########################################
    #find cluster median
    CL_median=lapply(CL,function(x){apply(expr[x,],2,median)})
    CL_median=data.frame(matrix(unlist(CL_median), ncol=length(CL_median[[1]]), byrow=T))
    CL_median=cbind(paste0("cluster",1:length(CL)),unique(CL_label),CL_median)
    names(CL_median)=c("cluster_id","label",antibodies)
    write.csv(CL_median,paste(outpath,"cluster_median.csv",sep="/"),row.names=F)
    cluster_summary=rbind(cluster_summary,data.frame("cluster_id"=CL_median$cluster_id,"label"=unique(CL_label)))

    ##### 5) record sample data#########################################
    SP_stat=NULL
    for(i in 1:length(CL)){
      for(j in unique(expr_sample)){
        w=intersect( which(expr_sample==j), CL[[i]] )
        frac=length(w)/nrow(expr)
        if(length(w)>1){M=apply(expr[w,],2,median)}
        if(length(w)==1){M=expr[w,]}
        if(length(w)==0){M=rep(NA,ncol(expr))}
        t1=c(paste0("cluster",i),CL_label[i],fcsFile,M,frac)
        SP_stat=rbind(SP_stat,t1)
      }
    }
    colnames(SP_stat)=c("cluster_id","label","fcs_files",antibodies,"fraction")
    write.csv(SP_stat,paste(outpath,"cluster_stats_in_each_sample.csv",sep="/"),row.names=F)

    ##### 6) plot the density plot############################
    if(ifPlot==T){
      filename=paste(outpath,"density_plot.pdf",sep="/")
      height=3*length(CL);width=3*length(antibodies)
      pdf(filename,width=width,height=height)
      par( mfcol = c(length(CL), length(antibodies) ) )
      for(i in 1:length(antibodies)){
        x_all=expr[,i]
        for(j in 1:length(CL)){
          b=seq(min(x_all),max(x_all), ((max(x_all)-min(x_all))/100) )
          hist(x_all,col=rgb(0, 0, 0, 0.2),xlab=antibodies[i],breaks=b,freq=T,border=F,main=paste0(antibodies[i]," of cluster ",j," : \n", CL_label[j]))
          hist(expr[CL[[j]],i],add=T,breaks=b,col=rgb(1, 0, 0, 0.5),freq=T,border=F)
          abline(v=cutoff[i])
          if(antibodies[i]%in%Ab_tri){abline(v=triS[,i])}
        }
      }
      dev.off()
    }

    ## 9) write.out cell id that belong to each clusters--------
    filename=paste(outpath,"cluster_cell_id.RData",sep="/")
    names(CL)=CL_label
    save(CL,file=filename)

}

flowHC=function(exprM,minimumClusterSizePercent=0.05){
  D=dist(exprM)
  HC=fastcluster::hclust(D,method="ward.D")
  CL=HC.assign(HC,minimumClusterSizePercent)
  return(CL)
}

skewness= function(v,x){
  c1=1*(x>v)
  x1=x[c1==0]
  x2=x[c1==1]
  
  t1=mean(x1)-median(x1)
  #t1=abs(t1)/(max(x1)-min(x1))*length(x1)
  t1=abs(t1)/sd(x1)*length(x1)
  
  t2=mean(x2)-median(x2)
  #t2=abs(t2)/(max(x2)-min(x2))*length(x2)
  t2=abs(t2)/sd(x2)*length(x2)
  return((t1+t2)/length(x))
}

flowSOM.MC=function(exprM,...){
  input=flowFrame(exprM)
  fSOM <- FlowSOM::ReadInput(input, transform = FALSE, scale = FALSE)
  fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse = 1:ncol(exprM), 
                            xdim = 10, ydim = 10)
  fSOM <- FlowSOM::BuildMST(fSOM)
  meta <- FlowSOM::metaClustering_consensus(fSOM$map$codes,k=40)
  CL <- meta[fSOM$map$mapping[, 1]]
  CL=lapply(unique(CL),function(x){which(CL==x)})
  return(CL)
}

