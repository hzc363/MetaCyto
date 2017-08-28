#' Find markers in a flow frame object
#'
#' A function that finds markers in a flow frame object.
#' @param fcsFrame A flow frame object.
#' @return Returns a vector of markers.
#' @details If the antibody name is available, the antibody name will be
#'   returned, otherwise, the channel name will be returned.
#' @examples
#' library(flowCore)
#' files=system.file("extdata","SDY420/ResultFiles/CyTOF_result",
#'                   package="MetaCyto")
#' files=list.files(files,pattern="fcs$",full.names=TRUE)[1]
#' fcs = read.FCS(files)
#' markers = markerFinder(fcs)
#' @importFrom flowCore pData parameters
#' @export
markerFinder=function(fcsFrame){
  Label=flowCore::pData(flowCore::parameters(fcsFrame));
  channels=Label$name
  antibodies=Label$desc
  antibodies=toupper(antibodies)
  channels=toupper(channels)
  antibodies=sapply(1:length(antibodies), function(i){
    if(is.na(antibodies[i])|antibodies[i]=="NA"){return(channels[i])}else{antibodies[i]}
  })
  return(antibodies)
}
