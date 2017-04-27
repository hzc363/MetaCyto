#' draw density plot for each cell cluster.
#'
#' A function that draws density plot for each cell cluster.
#' @param fcsFrame A flow Frame object returned from preprocessing function.
#' @param clusterList a list, each element should be a vector containing the IDs of all cells that belong to a cluster
#' @param cutoff a vector of cutoff values to bisect the distribution of each
#'   marker. The names of the vector should be the same as the marker names.
#' @param markerToPlot a vector specifying markers included in the plot. If NULL, all markers will be plotted.
#' @return NULL. The plot will show up automatically.
#' @details The plot can be very large, we suggest plotting it into a pdf jpeg file directly.
#' @examples
#' # Find fcs files
#' files=system.file("extdata","SDY420/ResultFiles/CyTOF_result",package="MetaCyto")
#' files=list.files(files,pattern="fcs$",full.names=TRUE)
#' # Preprocess
#' fcs = preprocessing(fcsFiles=files,assay ="CyTOF",b=1/8)
#' # Search clusters
#' cluster_list=searchCluster(fcsFrame=fcs,
#'                            clusterLabel=c("CD3+|CD8+","CD3-|CD19+"))
#' # plot density plot for clusters
#' densityPlot(fcs,
#'             clusterList=cluster_list$clusterList,
#'             cutoff=cluster_list$cutoff,
#'             markerToPlot=c("CD3","CD8","CD19"))
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
      hist(x_all,col=rgb(0, 0, 0, 0.5),xlab=markerToPlot[i],breaks=b,freq=TRUE,border=FALSE,main=paste0(markerToPlot[i]," of cluster ",j," : \n", CL_label[j]))
      hist(expr[CL[[j]],markerToPlot[i]],add=TRUE,breaks=b,col=rgb(1, 0, 0, 0.5),freq=TRUE,border=FALSE)
      if(markerToPlot[i]%in%names(cutoff)){abline(v=cutoff[markerToPlot[i]])}
    }
  }
}
