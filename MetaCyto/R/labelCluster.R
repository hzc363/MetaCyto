#' label each clusters as "+" or "-" or neutral for each marker
#'
#' A function that label each cluster as "+" or "-" or neutral for each marker
#' @param fcsFrame A flowFrame object.
#' @param clusterList a list, each element should be a vector containing the IDs of all cells that belongs to a cluster
#' @param excludeClusterParameters A vector specifiying the name of markers not to be used for labeling.
#' @param minPercent a number between 0 and 0.5. Used to specify the minimum percent of cells in positive and negative region after bisection. Keep it small to avoid bisecting uni-mode distributions.
#' @param labelQuantile a number between 0.5 and 1. Used to specify the minimum percent of a cluster required to be larger or smaller than the cutoff value for labeling.
#' @param cutoff a vector of cutoff values to bisect the distribution of each
#'   marker. The names of the vector should be the same as the marker names. If
#'   NULL, the cutoff value will be determined automatically.
#' @return Returns a list with two components: 1) clusterLabel, contains a vector of lables, each correspond to a cluster in clusterList. 2) cutoff, contains a vector of cutoff values used to bisect each marker.
#' @examples
#' # Find fcs files
#' files=system.file("extdata","SDY420/ResultFiles/CyTOF_result",package="MetaCyto")
#' files=list.files(files,pattern="fcs$",full.names=TRUE)
#' # Preprocess
#' fcs = preprocessing(fcsFiles=files,assay ="CyTOF",b=1/8)
#' # cluster using flowSOM.MC
#' cluster_list=flowSOM.MC(fcsFrame=fcs,excludeClusterParameters=c("Time","Cell_length"))
#' # label each clusters
#' cluster_label=labelCluster(fcsFrame=fcs,clusterList=cluster_list,
#'                            excludeClusterParameters=c("Time","Cell_length"))

#' @export
labelCluster= function(fcsFrame,
                       clusterList,
                       excludeClusterParameters=c("TIME"),
                       minPercent=0.05,
                       labelQuantile=0.95,
                       cutoff=NULL){

  #prepare exclude parameters
  excludeClusterParameters=toupper(excludeClusterParameters)
  excludeClusterParameters=union("SAMPLE_ID",excludeClusterParameters)

  # get marker names
  antibodies=markerFinder(fcsFrame)

  # Get expression matrix
  expr=flowCore::exprs(fcsFrame);
  colnames(expr)=antibodies

  # subset on columns
  w = (!antibodies%in%excludeClusterParameters) #& (!channels%in%excludeClusterParameters)
  antibodies=antibodies[w]
  expr=expr[,w,drop=FALSE]

  #find cutoff of each parameter
  if(is.null(cutoff)){
    cutoff=apply(X=expr,MARGIN=2,FUN=findCutoff)
  }

  #label each cluster
  CL_label=rep(NA,length(clusterList))
  no_label=sapply(1:length(antibodies),function(i){
    x1=quantile(expr[,i],minPercent)
    x2= quantile(expr[,i],(1-minPercent))
    if(cutoff[i]<x1|cutoff[i]>x2){return(TRUE)}else{return(FALSE)}
  })
  for(i in 1:length(clusterList)){
    l=c()
    for(j in 1:length(antibodies)){
      if(no_label[j]==TRUE){next}
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
