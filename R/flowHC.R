#' cluster cytometry data using hierarchical clustering
#'
#' A function that cluster cytometry data using hierarchical clustering.
#' @param fcsFrame a flow frame.
#' @param excludeClusterParameters A vector specifying the name of markers not to be used for clustering.
#' @param minimumClusterSizePercent a number between 0 and 1, used to specify the minimum size of a cluster relative to all events.
#' @return a list of clusters. Each cluster contains the ID of all cells that belong to the cluster.
#' @examples
#' # Find fcs files
#' files=system.file("extdata","SDY420/ResultFiles/CyTOF_result",package="MetaCyto")
#' files=list.files(files,pattern="fcs$",full.names=TRUE)
#' # Preprocess
#' fcs = preprocessing(fcsFiles=files,assay ="CyTOF",b=1/8)
#' # cluster using flowHC
#' cluster_list=flowHC(fcsFrame=fcs,excludeClusterParameters=c("Time","Cell_length"))
#' @export
#' @import fastcluster
flowHC=function(fcsFrame,excludeClusterParameters,minimumClusterSizePercent=0.05){
  antibodies=markerFinder(fcsFrame)
  excludeClusterParameters=toupper(excludeClusterParameters)
  excludeClusterParameters=union(excludeClusterParameters,c("TIME","SAMPLE_ID"))
  w=which(!antibodies%in%excludeClusterParameters)
  expr=flowCore::exprs(fcsFrame)
  expr=expr[,w,drop=FALSE]
  D=dist(expr)
  HC=fastcluster::hclust(D,method="ward.D")
  CL=HC.assign(HC,minimumClusterSizePercent)
  return(CL)
}
