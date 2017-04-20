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
                         b=1/150,
                         fileSampleSize=5000,
                         excludeTransformParameters=c("FSC-A","FSC-W","FSC-H","Time","Cell_length")){

  # loading libraries
  library(flowCore)
  library(dplyr)
  library(flowAI)

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
    fcs = read.flowSet(fcs_files, transformation = "linearize", alter.names = FALSE)
    if(!is.null(fileSampleSize)){
      fcs = fsApply(fcs, function(f){
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
      version = fsApply(fcs, function(frame) {
        return(keyword(frame, "FCSversion")[[1]])
      })
      unique_version = as.numeric(unique(as.vector(version)))

      # compensation, transformation
      ## if the file version is 2.0, then code will log transform the data
      if (unique_version == 2) {
        trans = logTransform()
        translist = transformList(biomarker_vector,trans)
        fcs = flowCore::transform(fcs, translist)
      } else if (unique_version == 3) {

        ### function first checks if the spill matrix in the flowframe is an actual spill matrix or an identity matrix
        if (is.null(keyword(fcs[[1]], "SPILL")[[1]]) == FALSE) {
          check = fsApply(fcs, function(x) {
            result = isSymmetric(keyword(x, "SPILL")[[1]])
          })
          ### if all flowframes have spill matrices then my function
          ### will apply them to each flowframe in the flowset for compensation
          if (unique(check)[[1]] == FALSE) {
            fcs = fsApply(fcs, function(x) {
              new_frame = compensate(x, keyword(x, "SPILL")[[1]])
              return(new_frame)
            })
          }
        }


        trans = arcsinhTransform(transformationId = "defaultArcsinhTransform", a = 0, b = b, c = 0)
        translist = transformList(biomarker_vector, trans)
        fcs = flowCore::transform(fcs, translist)
      }
    } else if (assay == "CyTOF") {
      w=which(!grepl(excludeTransformParameters,flowCore::colnames(fcs),ignore.case = T))
      biomarker_vector = flowCore::colnames(fcs)[w]

      trans = arcsinhTransform(transformationId = "defaultArcsinhTransform", a = 0, b = b, c = 0)
      translist = transformList(biomarker_vector, trans)
      fcs = flowCore::transform(fcs, translist)
    }

    # # 4) Remove inconsistency from parameter names (convert to uppercase and remove whitespace)
    # #Define function
    # f = function(frame) {
    #   # get channel names
    #   channels = as.vector(pData(parameters(frame))[,1])
    #   # get antibody names
    #   antibodies = as.vector(pData(parameters(frame))[,2])
    #   # combine both into a vector
    #   parameters = c(antibodies,channels)
    #   # remove NA values
    #   parameters = parameters[!is.na(parameters)]
    # 
    #   #change description values
    #   # iterating over channel and antibody names in parameters vector
    #   for (i in parameters) {
    #     # within the description slot, function locates where antibody/channel is located and returns the index
    #     location = which(description(frame)==i)
    #     # for each location function will convert name to uppercase and remove whitespace
    #     description(frame)[[location]] = toupper(description(frame)[[location]])
    #     description(frame)[[location]] = gsub('\\s+', '',description(frame)[[location]])
    #   }
    # 
    #   # within the desription slot, change spill matrix column names
    #   if(!is.null(description(frame)[["SPILL"]])){
    #     colnames(description(frame)[["SPILL"]]) = toupper(colnames(description(frame)[["SPILL"]]))
    #     colnames(description(frame)[["SPILL"]]) = gsub('\\s+', '',colnames(description(frame)[["SPILL"]]))
    #   }
    # 
    #   # within the exprs slot, change exprs column names
    #   flowCore::colnames(frame) = toupper(flowCore::colnames(frame))
    #   flowCore::colnames(frame) = gsub('\\s+', '',flowCore::colnames(frame))
    # 
    #   # within the parameter slot, change parameter names
    #   param = pData(parameters(frame))
    #   param$desc = toupper(param$desc)
    #   param$desc = gsub('\\s+', '',param$desc)
    #   pData(parameters(frame)) = param
    #   return(frame)
    # }
    # # iterating through each flowframe
    # fcs = fsApply(fcs, f)
    
    # Quality check
    # fcs=flow_auto_qc(fcs0, remove_from = "FM", output = 1, timeCh = NULL,
    #                 second_fractionFR = 0.1, alphaFR = 0.01, decompFR = TRUE,
    #                 ChRemoveFS = c("FSC", "SSC"), outlierFS = FALSE, pen_valueFS = 200,
    #                 max_cptFS = 3, ChFM = NULL, sideFM = "both", neg_valuesFM = 1,
    #                 html_report = FALSE, mini_report = F, fcs_QC = F,
    #                 fcs_highQ = FALSE, fcs_lowQ = FALSE, folder_results = "resultsQC")
    

    # 5) outputting resuts
    write.flowSet(fcs, outdir = paste0(outpath,'/',std))

    #a function to extract parameter names
    fp = function(frame) {
      # get channel names
      channels=as.vector(pData(parameters(frame))[,1])
      # get antibody names
      antibodies=as.vector(pData(parameters(frame))[,2])
      antibodies=sapply(1:length(antibodies), function(i){
        if(is.na(antibodies[i])|antibodies[i]=="NA"){return(channels[i])}else{antibodies[i]}
      })

      channels = paste(channels,collapse='|' )
      antibodies = paste(antibodies,collapse='|' )
      return(c(channels,antibodies))
    }
    param= fsApply(fcs, fp)
    param = toupper(param)
    fcs_param=rbind(fcs_param,data.frame("fcs_names"=rownames(param),"channels"=param[,1],"antibodies"=param[,2]))
  }#end of each std
  inputMeta=left_join(fcs_param,inputMeta,by="fcs_names")
  #inputMeta=select(inputMeta,-fcs_files)
  write.csv(inputMeta,paste0(outpath,"/processed_sample_summary.csv"),row.names=F)
  cat("Preprocess result stored in the folder:",outpath, "\n")
}#end of function

autoCluster= function(preprocessOutputFolder,
                      outpath="cluster_output",
                      ###tuning parameters##
                      excludeClusterParameters=c("TIME"),
                      minimumClusterSizePercent=0.05,
                      fileSampleSize=1000,
                      labelQuantile=0.8){

  ## load libraries
  library(flowCore)
  library(dplyr)
  library(fastcluster)


  #read the output from preprocessing
  inputMeta=read.csv(paste0(preprocessOutputFolder,'/processed_sample_summary.csv'),stringsAsFactors=F)
  #create output foler
  dir.create(file.path(outpath),recursive=T)

  #prepare exclude parameters
  excludeClusterParameters=toupper(excludeClusterParameters)
  cluster_summary=NULL
  for(std in unique(inputMeta$study_id)){
    cat("Clustering , study ID = ",std, "\n")

    ##### 1) read sample files for each study
    fcs_files=subset(inputMeta$fcs_names,inputMeta$study_id==std)
    fcs_files=paste(preprocessOutputFolder,std,fcs_files,sep="/")
    fcs=read.flowSet(fcs_files, transformation = "linearize", alter.names = FALSE)
    # make sure the fcs file antibody names are the same as the preprocessed output
    antibodies=subset(inputMeta$antibodies,inputMeta$study_id==std)[1]
    antibodies=strsplit(antibodies,"\\|")[[1]]

    ##### 2) subset the cells in fcs
    # Get expression matrix
    expr=fsApply(fcs,function(x){
      v=exprs(x);
      v=v[sample(1:nrow(v),fileSampleSize,replace=T),]
    })
    colnames(expr)=antibodies
    expr_sample=ceiling(seq_along(1:nrow(expr))/fileSampleSize)


    # subset on columns
    w=!antibodies%in%excludeClusterParameters
    antibodies=antibodies[w]
    if(length(antibodies)<2){next}
    expr=expr[,w]
    expr_scale=scale(expr,center=F,scale=T)

    # 3) cluster the subset samples, find clusters large enough
    ############################################################
    max_cell_N=20000
    if(nrow(expr_scale)>max_cell_N){
      ind=sample(1:nrow(expr_scale),max_cell_N,replace=F)

      D=dist(expr_scale[ind,])
    }else{
      D=dist(expr_scale)
    }
    HC=fastcluster::hclust(D,method="ward.D")
    CL=HC.assign(HC,minimumClusterSizePercent)

    # 4) label the clusters:
    ############################################################
    #find cutoff of each parameter
    cutoff=apply(expr,2, findCutoff)

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

    ##### 5) merge clusters
    #######################################################
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

    # 6) assign cluster to all cells
    #########################################
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

    ### 6) record cluster data
    #########################################
    dir.create(file.path(paste(outpath,std,sep="/")))
    #find cluster median
    CL_median=lapply(CL,function(x){apply(expr[x,],2,median)})
    CL_median=data.frame(matrix(unlist(CL_median), ncol=length(CL_median[[1]]), byrow=T))
    CL_median=cbind(std,paste0("cluster",1:length(CL)),unique(CL_label),CL_median)
    names(CL_median)=c("study","cluster_id","label",antibodies)
    write.csv(CL_median,paste(outpath,std,"cluster_median.csv",sep="/"),row.names=F)
    cluster_summary=rbind(cluster_summary,data.frame("Study"=std,"cluster_id"=CL_median$cluster_id,"label"=unique(CL_label)))

    ### 7) record sample data
    #########################################
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


    ## 8) plot the density plot
    ############################
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
  write.csv(cluster_summary,paste(outpath,"cluster_list.csv",sep="/"),row.names=F)
}

searchCluster=function(preprocessOutputFolder,
                       outpath="search_output",
                       clusterLabel,
                       fileSampleSize=1000,
                       QCpercent=0.0025){

  ## load libraries
  library(flowCore)
  library(dplyr)

  #read the output from preprocessing
  inputMeta=read.csv(paste0(preprocessOutputFolder,'/processed_sample_summary.csv'),stringsAsFactors=F)
  #create output foler
  dir.create(file.path(outpath),recursive=T)

  #prepare exclude parameters
  clusterLabel=toupper(clusterLabel)
  cluster_summary=NULL
  for(std in unique(inputMeta$study_id)){
    cat("Searching , study ID = ",std, "\n")

    ##### 1) read sample files for each study
    ############################################################
    fcs_files=subset(inputMeta$fcs_names,inputMeta$study_id==std)
    fcs_files=paste(preprocessOutputFolder,std,fcs_files,sep="/")
    fcs=read.flowSet(fcs_files, transformation = "linearize", alter.names = FALSE)
    # make sure the fcs file antibody names are the same as the preprocessed output
    antibodies=subset(inputMeta$antibodies,inputMeta$study_id==std)[1]
    antibodies=strsplit(antibodies,"\\|")[[1]]

    ##### 2) subset the cells in fcs
    ############################################################
    # Get expression matrix
    expr=fsApply(fcs,function(x){
      v=exprs(x);
      v=v[sample(1:nrow(v),fileSampleSize,replace=T),]
    })
    row_RM=lapply(as.data.frame(expr),function(x){which(x<quantile(x,QCpercent) |x>quantile(x,(1-QCpercent)) )})
    row_RM=Reduce(union,row_RM)
    if(length(row_RM)>0){expr=expr[-row_RM,]}
    colnames(expr)=antibodies
    expr_sample=ceiling(seq_along(1:nrow(expr))/fileSampleSize)

    # 3) bisect each marker
    ############################################################
    #find cutoff of each parameter
    cutoff=apply(expr,2, findCutoff)

    #find cluster whose markers are included in the study
    CL_label=clusterLabel
    CL_label2=strsplit(clusterLabel,split="&|\\|")
    in_study=sapply(CL_label2,function(x){
      Ab_in_label=gsub("\\+$|-$","",x)
      if(length(setdiff(Ab_in_label,antibodies))>0){return(F)}else{return(T)}
    })
    CL_label=CL_label[in_study]
    CL_label2=CL_label2[in_study]
    if(length(CL_label)<1){next}

    #find each cluster
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
    CL_count=sapply(CL,length)
    CL=CL[CL_count>1]
    CL_label=CL_label[CL_count>1]
    if(length(CL)<1){next}

    ### 6) record cluster data
    #########################################
    dir.create(file.path(paste(outpath,std,sep="/")))
    #find cluster median
    CL_median=lapply(CL,function(x){apply(expr[x,],2,median)})
    CL_median=data.frame(matrix(unlist(CL_median), ncol=length(CL_median[[1]]), byrow=T))
    CL_median=cbind(std,paste0("cluster",1:length(CL)),unique(CL_label),CL_median)
    names(CL_median)=c("study","cluster_id","label",antibodies)
    write.csv(CL_median,paste(outpath,std,"cluster_median.csv",sep="/"),row.names=F)
    cluster_summary=rbind(cluster_summary,data.frame("Study"=std,"cluster_id"=CL_median$cluster_id,"label"=unique(CL_label)))

    ### 7) record sample data
    #########################################
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


    ## 8) plot the density plot
    ############################
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
  write.csv(cluster_summary,paste(outpath,"cluster_list.csv",sep="/"),row.names=F)
}

fcsInfoParser =function(metaData,
                        studyFolder,
                        fcsCol="ZBXFN",
                        assay=c("FCM", "CyTOF")){
  library(flowCore)
  cat("Parsing the meta-data from ImmPort\n")
  metaData=unique(metaData[,fcsCol,drop=FALSE])


  if(assay=="FCM"){
    fcs_files=paste(studyFolder,"ResultFiles/Flow_cytometry_result",metaData[,fcsCol],sep="/")
  }else{
    fcs_files=paste(studyFolder,"ResultFiles/CyTOF_result",metaData[,fcsCol],sep="/")
  }
  fsc_markers=sapply(fcs_files,function(x){
    tryCatch({
      fcs=read.FCS(x)
      panel=paste( pData(parameters(fcs))$name, pData(parameters(fcs))$desc, sep="-")
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


collectData=function(files,longform=T){
  library(tidyr)
  all_data=NULL
  nms=NULL
  for(fn in files){
    #read output data from "autoCluster" function
    cluster_stat=read.csv(fn,stringsAsFactors=F,check.names=F)
    if(is.null(nms)){nms=colnames(cluster_stat)}
    w=which(nms=="fcs_files")+1
    if(longform==T){
      cluster_stat=gather(cluster_stat,parameter_name, value ,w:ncol(cluster_stat))
    }else(colnames(cluster_stat)=nms)
    all_data=rbind(all_data,cluster_stat)
  }
  return(all_data)
}

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

panelSummary=function(panelInfo,folder,cluster=T,plotImage=T){
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
  ab_table=spread(ab_list,key=antibodies,value=value,fill=0)
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
    pdf(paste0(folder,"/panel_summary.pdf"),width=ncol(ab_table)*3,height=nrow(ab_table)*1)
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

metaAnalysis=function(value,variableOfInterst,otherVariables,
                      studyID,data,CILevel,main,
                      ifScale=c(T,F)){
  library(metafor)
  study_result=NULL
  for(std in unique(data[,studyID])){
    sub_data=subset(data,data[,studyID]==std)
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
                  xlab="Effect Size", mlab="RE Model for All Studies")
  t1=data.frame("Summary",res$b[1],res$se,res$zval,res$pval,res$ci.lb,res$ci.ub,
                sum(study_result$N))
  names(t1)=names(study_result)
  study_result=rbind(study_result,t1)
  rownames(study_result)=NULL
  return(study_result)
}

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

nameUpdator=function(oldNames,newNames,files){
  for(fn in files){
    TB=read.csv(fn,stringsAsFactors=F,check.names=F)
    for(i in 1:length(newNames)){
      TB$antibodies=sapply(TB$antibodies,function(x){
        gsub(old_names[i],new_names[i],x,fixed=T)
      })
    }
    write.csv(TB,fn,row.names=F)
  }
}


# Define back-end functions----------
findCutoff=function(x){
  library(cluster)
  if(length(x)>2000){x=sample(x,2000)}
  x=x[x<quantile(x,0.98)&x>quantile(x,0.02)]
  x_ori=x
  x=spaceRemove(x)
  valley=sapply(seq(0.01,0.99,length.out=100),function(q){quantile(x,q)})
  D=dist(x)
  sil=sapply(valley,function(v){
    cluster=1*(x>v)+1
    if(length(unique(cluster))<2){return(-2)}
    ss <- silhouette(cluster,D)
    return(mean(ss[, 3]))
  })
  #valley=valley[which.max(sil)]
  valley=quantile(x_ori,(which.max(sil))/100)
  # km=kmeans(x,2)
  # t1=c(range(x[km$cluster==1]),range(x[km$cluster==2]))
  # t1=sort(t1)
  # valley=mean(t1[2:3])
  return(valley)
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
  for(k in 2:nrow(HC$merge)){
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

spaceRemove=function(x,N.break=length(x)/20){
  GRID=seq(min(x),max(x),length.out=N.break)
  #GRID=GRID[-c(1,length(GRID))]
  D_GRID=diff(GRID)
  X_bin=.bincode(x,GRID,include.lowest=T)
  N_each=rep(0,length(D_GRID))
  N_each[sort(unique(X_bin))]=table(sort(X_bin))
  
  D_GRID_new=D_GRID*(N_each>0)
  GRID_new=GRID[1]
  for(i in 1:length(D_GRID)){
    GRID_new[i+1]=GRID_new[i]+D_GRID_new[i]
  }
  
  x_new=lapply(unique(X_bin),function(i){ x[X_bin==i]-GRID[i]+GRID_new[i] })
  x_new=unlist(x_new)
  return(x_new)
}





