#' Collect and combine data from multiple csv files of the same format
#'
#' A function that collect and combine data from multiple csv files of the same
#' format.
#' @param files A vector containing the paths of csv files to be combined.
#' @param longform True or False. Used to specify if the table in each csv file
#'   should be converted into long form before combining.
#' @return A dataframe containing conbined information from multiple csv files.
#' @examples
#' # find all the files we want to combine
#' fn=system.file("extdata","",package="MetaCyto")
#' fn=list.files(fn,pattern="cluster_stats_in_each_sample",full.names=TRUE)
#' # Comine the data
#' all_data = collectData(fn,longform=TRUE)
#' @importFrom tidyr gather
#' @export
collectData=function(files,longform=TRUE){
  all_data=NULL
  nms=NULL
  for(fn in files){
    #read output data from "autoCluster" function
    cluster_stat=read.csv(fn,stringsAsFactors=FALSE,check.names=FALSE)
    if(is.null(nms)){nms=colnames(cluster_stat)}
    w=which(nms=="fcs_names")+1
    if(longform==TRUE){
      cluster_stat=tidyr::gather(cluster_stat,parameter_name, value ,w:ncol(cluster_stat))
    }else(colnames(cluster_stat)=nms)
    all_data=rbind(all_data,cluster_stat)
  }
  return(all_data)
}
