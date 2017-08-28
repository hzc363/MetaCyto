#' Find cutoff for a 1D distribution
#'
#' A function that finds cutoff for a 1D distribution.
#' @param x A vector of values.
#' @param returnSil Logic, used to specify if the max average silhouette is
#'   returned
#' @param useBL Logic, used to specify if outliers should be ignored
#' @param minX A numerical value, used to specify the min value allowed for the
#'   cutoff.
#' @return If returnSil=F, returns a single cutoff value. Otherwise, returns a
#'   list containing the cutoff value and the max average silhouette
#' @examples
#' x=c(rnorm(1000),rnorm(1000,5))
#' findCutoff(x)
#' @importFrom cluster silhouette
#' @export
findCutoff=function(x,returnSil=FALSE,useBL=TRUE,minX=0){
  if(length(x)>2000){x=sample(x,2000)}
  if(useBL==TRUE){x=baselineCut(x)}
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
  if(returnSil==FALSE){return(valley)}else{return(c("cuoff"=valley,"sil"=max(sil)))}
}
