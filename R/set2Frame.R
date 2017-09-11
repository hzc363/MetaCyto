#' Combine cells in a flow set into a flow frame.
#'
#' A function that combines cells in a flow set into a flow frame.
#' @param flowSet A flow set object
#' @return Returns a flowFrame object. All cells from flow set are combined into
#'   one flow frame. A new parameter, sample_id, is introduced to indicate the
#'   origin of each cell.
#' @examples
#' library(flowCore)
#' files=system.file("extdata","SDY420/ResultFiles/CyTOF_result",
#'                   package="MetaCyto")
#' files=list.files(files,pattern="fcs$",full.names=TRUE)
#' flow_set = read.flowSet(files)
#' flow_frame = set2Frame(flow_set)
#' @importFrom flowCore fsApply exprs pData flowFrame parameters
#' @export
set2Frame=function(flowSet){
  expr=flowCore::fsApply(flowSet,function(x){
    v=flowCore::exprs(x);
    return(v)
  })
  eventN=flowCore::fsApply(flowSet,function(x){
    v=flowCore::exprs(x);
    n=nrow(v)
    return(n)
  })
  Label=flowCore::fsApply(flowSet,function(x){
    v=flowCore::pData(flowCore::parameters(x));
    return(v)
  })
  #annotate expr colnames
  channels=Label[[1]]$name
  antibodies=Label[[1]]$desc
  antibodies=toupper(antibodies)
  channels=toupper(channels)
  antibodies=sapply(1:length(antibodies), function(i){
    if(is.na(antibodies[i])|antibodies[i]=="NA"){return(channels[i])}else{antibodies[i]}
  })
  colnames(expr)=antibodies
  #add sample_id
  sample_id=lapply(1:length(eventN),function(i){rep(i,eventN[i])})
  sample_id=unlist(sample_id)
  expr=cbind(expr,"sample_id"=sample_id)
  #return flowFrame
  fFrame=flowCore::flowFrame(expr)
  rownames(fFrame@parameters@data) = paste0("$P",1:nrow(fFrame@parameters@data))
  return(fFrame)
}
