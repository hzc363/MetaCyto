#' Collect sample information for fcs files in a study from ImmPort.
#'
#' A function that collects sample information for fcs files in a study from
#' ImmPort.
#' @param metaData A data frame. Must contain a column listing the names of fcs
#'   files included in the study.
#' @param studyFolder Path of the directory containing all the files of a study
#'   from ImmPort.
#' @param fcsCol A string specifying the name of the column in metaData that
#'   lists fcs files included in the study.
#' @param assay Either "FCM" or "CyTOF" to indicate the type of cytometry data.
#' @param attrCol A vector of column names. Used to specify the information
#'   about each cytometry the user wish to include in the analysis.
#' @return A dataframe containing sample information.
#' @examples
#' fn=system.file("extdata",
#'                "SDY736/SDY736-DR19_Subject_2_Flow_cytometry_result.txt",
#'                package="MetaCyto")
#' meta_data=read.table(fn,sep='\t',header=TRUE,check.names=FALSE)
#' # Find the AGE, GENDER info from selected_data
#' fn=system.file("extdata","SDY736",package="MetaCyto")
#' sample_info_SDY736=sampleInfoParser(metaData=meta_data,
#'                                     studyFolder=fn,
#'                                     assay="FCM",
#'                                     fcsCol="File Name",
#'                                     attrCol=c("Subject Age","Gender"))
#' @export
sampleInfoParser =function(metaData,
                           studyFolder,
                           fcsCol="ZBXFN",
                           assay="FCM",
                           attrCol){
  if(assay=="FCM"){
    fcs_files=paste(studyFolder,"ResultFiles/Flow_cytometry_result",metaData[,fcsCol],sep="/")
  }else{
    fcs_files=paste(studyFolder,"ResultFiles/CyTOF_result",metaData[,fcsCol],sep="/")
  }
  inputMeta=data.frame("fcs_files"=fcs_files)
  inputMeta=cbind(inputMeta,metaData[,attrCol])
  inputMeta=unique(inputMeta)
  inputMeta$fcs_files=as.character(inputMeta$fcs_files)
  return(inputMeta)
}
