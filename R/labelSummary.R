#' Make a summary for the labels (cell populations) identified in different
#' cytometry panels.
#'
#' A function that summarizes the labels (cell populations) identified in
#' different cytometry panels.
#' @param allData A table containing the summary statistics of cell populations.
#'   Often is the output of the function "collectData".
#' @param minStudy Minimum number of cytometry panels a label should appear in.
#'   Set it >1 to find cell populations identified in more than 1 cytometry
#'   panel for meta-analysis.
#' @return A data frame summarizing the labels identified in different cytometry
#'   panels.
#' @examples
#' fn=system.file("extdata","",package="MetaCyto")
#' files=list.files(fn,pattern="cluster_stats_in_each_sample",
#'                  recursive=TRUE,full.names=TRUE)
#' fcs_stats=collectData(files,longform=TRUE)
#' label_summary = labelSummary(allData=fcs_stats,minStudy=2)
#' @export
labelSummary=function(allData,minStudy=2){
  result=NULL
  for(label in unique(allData$label)){
    mid_data=allData[allData$label==label,]
    for(param in unique(allData$parameter_name)){
      sub_data=mid_data[mid_data$parameter_name==param,]
      if(nrow(sub_data)<1){next}
      N_study=length(unique(sub_data$study))
      if(N_study<minStudy){next}
      t1=data.frame("label"=label,"parameter_name"= param,"number_of_study"=N_study,"sample_size"=nrow(sub_data),
                    "studies"=paste(unique(sub_data$study),collapse="|"), row.names=NULL )
      result=rbind(result,t1)
    }
  }
  return(result)
}
