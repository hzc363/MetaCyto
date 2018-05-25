#' Cluster cytometry data using FlowSOM
#'
#' A function that cluster cytometry data using FlowSOM.
#' @param fcsFrame A flow frame.
#' @param excludeClusterParameters A vector specifying the name of markers not
#'   to be used for clustering.
#' @param xdim An integer, defines the width of SOM
#' @param ydim An integer, defines the height of SOM
#' @param k An integer, defines the number of clusters to be identified by
#'   flowSOM.
#' @param seed Set seed for reproducible results.
#' @return A list of clusters. Each cluster contains the ID of all cells that
#'   belong to the cluster.
#' @examples
#' # Find fcs files
#' files=system.file("extdata",
#'                   "SDY420/ResultFiles/CyTOF_result",package="MetaCyto")
#' files=list.files(files,pattern="fcs$",full.names=TRUE)
#' # Preprocess
#' fcs = preprocessing(fcsFiles=files,assay ="CyTOF",b=1/8)
#' # cluster using flowSOM.MC
#' cluster_list=flowSOM.MC(fcsFrame=fcs,
#'                         excludeClusterParameters=c("Time","Cell_length"))
#' @importFrom FlowSOM ReadInput BuildSOM BuildMST metaClustering_consensus
#' @export
flowSOM.MC=function(fcsFrame,excludeClusterParameters,
                    xdim=10,ydim=10,k=40,seed = NULL){
  set.seed(seed)
  antibodies=markerFinder(fcsFrame)
  excludeClusterParameters=toupper(excludeClusterParameters)
  excludeClusterParameters=union(excludeClusterParameters,c("TIME","SAMPLE_ID"))
  w=which(!antibodies%in%excludeClusterParameters)
  fSOM <- FlowSOM::ReadInput(fcsFrame, transform = FALSE, scale = FALSE)
  fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse = w,
                            xdim = xdim, ydim = ydim)
  fSOM <- FlowSOM::BuildMST(fSOM)
  meta <- FlowSOM::metaClustering_consensus(fSOM$map$codes,k=k,seed = seed)
  CL <- meta[fSOM$map$mapping[, 1]]
  CL=lapply(unique(CL),function(x){which(CL==x)})
  return(CL)
}
