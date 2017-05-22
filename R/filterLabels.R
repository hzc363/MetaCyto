#' Filter cluster labels
#'
#' A function that filter cluster labels.
#' @param labels A vector containing labels for cell clusters
#' @param minPlus An integer, used to specify the minimum  number of "+" a label should contain.
#' @param minMarker An integer, used to specify the minimum  number of markers a label should contain.
#' @param maxMarker An integer, used to specify the max number of markers a label should contain.
#' @return Returns a vector of labels that pass through the filter.
#' @examples
#' labels= c("CD3-|CD4-|CD8-","CD3+|CD4+|CD8-","CD3+|CD4-|CD8+","CD3+|CD4-|CD8+|CCR7+|CD45RA-|CCR6-")
#' labels=filterLabels(labels=labels,minPlus=1,minMarker=2,maxMarker=5)
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
