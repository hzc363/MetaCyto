#' Search for clusters using pre-defined labels in cytometry data from different
#' studies in batch
#'
#' A function that searches for clusters using pre-defined labels in cytometry
#' data from different studies in batch.
#' @param preprocessOutputFolder Directory where the pre-processed results are
#'   stored.
#' @param outpath A string indicating the directory the results should be
#'   written to.
#' @param clusterLabel A vector containing labels, such as c("CD3+|CD4+|CD8-")
#' @param ifPlot True or False. Used to specify if a the density plot for each
#'   cluster should be plotted
#' @return Results will be written to the outpath folder
#' @details The function writes out the summary statistics for each cluster. A
#'   separate directory will be created for each study.
#' @examples
#' #get meta-data
#' fn=system.file("extdata","fcs_info.csv",package="MetaCyto")
#' fcs_info=read.csv(fn,stringsAsFactors=FALSE,check.names=FALSE)
#' fcs_info$fcs_files=system.file("extdata",fcs_info$fcs_files,
#'                                package="MetaCyto")
#' # Make sure the transformation parameter "b" and the "assay" argument
#' # are correct of FCM and CyTOF files
#' b=assay=rep(NA,nrow(fcs_info))
#' b[grepl("CyTOF",fcs_info$study_id)]=1/8
#' b[grepl("FCM",fcs_info$study_id)]=1/150
#' assay[grepl("CyTOF",fcs_info$study_id)]="CyTOF"
#' assay[grepl("FCM",fcs_info$study_id)]="FCM"
#' # preprocessing
#' preprocessing.batch(inputMeta=fcs_info,
#'                     assay=assay,
#'                     b=b,
#'                     outpath="Example_Result/preprocess_output",
#'                     excludeTransformParameters=c("FSC-A","FSC-W","FSC-H",
#'                     "Time","Cell_length"))
#' # Make sure marker names are consistant in different studies
#' files=list.files("Example_Result",pattern="processed_sample",recursive=TRUE,
#'                  full.names=TRUE)
#' nameUpdator("CD8B","CD8",files)
#' # find the clusters
#' cluster_label=c("CD3-|CD19+","CD3+|CD4-|CD8+")
#' searchCluster.batch(preprocessOutputFolder="Example_Result/preprocess_output",
#'                     outpath="Example_Result/search_output",
#'                     clusterLabel=cluster_label)
#' @importFrom flowCore read.FCS exprs flowFrame
#' @export
searchCluster.batch=function(preprocessOutputFolder,
                             outpath="search_output",
                             clusterLabel,
                             ifPlot=TRUE){

  #read the output from preprocessing
  inputMeta=read.csv(file.path(preprocessOutputFolder,'processed_sample_summary.csv'),stringsAsFactors=FALSE)
  #create output foler
  dir.create(file.path(outpath),recursive=TRUE)

  #prepare exclude parameters
  clusterLabel=toupper(clusterLabel)
  cluster_summary=NULL
  for(std in unique(inputMeta$study_id)){
    cat("Searching , study ID = ",std, "\n")
    dir.create(file.path(outpath,std))
    ##### 1) read sample files for each study############################################################
    fcs_files=file.path(preprocessOutputFolder,paste0(std,".fcs"))
    fcs=flowCore::read.FCS(fcs_files,truncate_max_range=FALSE)

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
    write.csv(CL_stats,file.path(outpath,std,"cluster_stats_in_each_sample.csv"),row.names=FALSE)

    ##### 6) plot the density plot############################
    if(ifPlot==TRUE){
      filename=file.path(outpath,std,"density_plot.pdf")
      height=3*length(CL);width=3*length(antibodies)
      pdf(filename,width=width,height=height)
      densityPlot(fcs,CL,cutoff=sC$cutoff,markerToPlot = names(sC$cutoff))
      dev.off()
    }
  }#end of each study
}
