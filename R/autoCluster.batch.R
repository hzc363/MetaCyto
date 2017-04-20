#' Cluster the preprocessed fcs files from different studies in batch
#'
#' A function that clusters the pre-processed fcs files from different studies in batch.
#' @param preprocessOutputFolder Directory where the preprocessed results are stored. Should be the same with the outpath argument in preprocessing.batch function.
#' @param excludeClusterParameters A vector specifiying the name of markers not to be used for clustering and labeling. Typical example includes: Time, cell_length.
#' @param labelQuantile A number between 0.5 and 1. Used to specify the minimum percent of cells in a cluster required to be larger or smaller than the cutoff value for labeling.
#' @param clusterFunction The name of unsuperviesed clustering function the user wish to use for clustering the cells. The default is "flowSOM.MC". The first argument of the function must take a flow frame, the second argument of the function must take a vector of excludeClusterParameters. The function must returns a list of clusters containing cell IDs. flowSOM.MC and flowHC are implemented in the package. For other methods, please make your own wrapper functions.
#' @param minPercent a number between 0 and 0.5. Used to specify the minimum percent of cells in positive and negative region after bisection. Keep it small to avoid bisecting uni-mode distributions.
#' @param ... pass arguments to labelCluster and clusterFunction
#' @return a vector of labels identified in the cytometry data.
#' @examples
#' #get meta-data
#' fn=system.file("extdata","fcs_info.csv",package="MetaCyto")
#' fcs_info=read.csv(fn,stringsAsFactors=FALSE,check.names=FALSE)
#' fcs_info$fcs_files=system.file("extdata",fcs_info$fcs_files,package="MetaCyto")
#' # make sure the transformation parameter "b" and the "assay" argument are correct of FCM and CyTOF files
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
#'                     excludeTransformParameters=c("FSC-A","FSC-W","FSC-H","Time","Cell_length"))
#' # Make sure marker names are consistant in different studies
#' files=list.files("Example_Result",pattern="processed_sample",recursive=TRUE,full.names=TRUE)
#' nameUpdator("CD8B","CD8",files)
#' # find the clusters
#' excludeClusterParameters=c("FSC-A","FSC-W","FSC-H","SSC-A","SSC-W","SSC-H","Time",
#'                           "CELL_LENGTH","DEAD","DNA1","DNA2")
#' cluster_label=autoCluster.batch(preprocessOutputFolder="Example_Result/preprocess_output",
#'                                 excludeClusterParameters=excludeClusterParameters,
#'                                 labelQuantile=0.95,
#'                                 clusterFunction=flowHC)
#' @export
autoCluster.batch= function(preprocessOutputFolder,
                            excludeClusterParameters=c("TIME"),
                            labelQuantile=0.9,
                            clusterFunction=flowSOM.MC,
                            minPercent=0.05,...){
  #read the output from preprocessing
  inputMeta=read.csv(paste0(preprocessOutputFolder,'/processed_sample_summary.csv'),stringsAsFactors=FALSE)
  #create output foler

  #prepare exclude parameters
  excludeClusterParameters=toupper(excludeClusterParameters)
  all_labels=NULL
  for(std in unique(inputMeta$study_id)){
    cat("Clustering , study ID = ",std, "\n")

    ##### 1) read sample files for each study############################################################
    fcs_files=paste0(preprocessOutputFolder,"/",std,".fcs")
    fcs=flowCore::read.FCS(fcs_files,truncate_max_range=FALSE)

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
    expr=expr[,w,drop=FALSE]
    expr_scale=scale(expr,center=FALSE,scale=TRUE)
    fcs=flowCore::flowFrame(expr_scale)

    CL=clusterFunction(fcs,excludeClusterParameters,...)
    CL_label=labelCluster(fcs,CL,excludeClusterParameters,
                          labelQuantile=labelQuantile,
                          minPercent=minPercent,...)
    all_labels=union(all_labels,CL_label$clusterLabel)
  }#end of each study
  return(all_labels)
}

