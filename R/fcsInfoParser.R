#' organize fcs files in a study from ImmPort into panels
#'
#' A function that organizes fcs files in a study from ImmPort into panels.
#' @param metaData a data frame. Must contain a column listing the names of fcs
#'   files included in the study.
#' @param studyFolder Path of the directory containing all the files of a study
#'   from ImmPort.
#' @param fcsCol a string specifying the name of the column in metaData that
#'   lists fcs files included in the study.
#' @param assay Either "FCM" or "CyTOF" to indicate the type of cytometry data.
#' @return A dataframe containing 2 columns: a column called "fcs_files" that
#'   contains the location (relative to the working directory) of each fcs file
#'   on the hard drive and a column called "study_id" that specify what study
#'   each fcs file belongs to.
#' @examples
#' fn=system.file("extdata","SDY736/SDY736-DR19_Subject_2_Flow_cytometry_result.txt",package="MetaCyto")
#' meta_data=read.table(fn,sep='\t',header=TRUE)
#' # Organize fcs file into panels
#' fn=system.file("extdata","SDY736",package="MetaCyto")
#' fcs_info_SDY736=fcsInfoParser(metaData=meta_data,
#'                               studyFolder=fn,
#'                               fcsCol="FILE_NAME",
#'                               assay="FCM")
#' @export
fcsInfoParser =function(metaData,
                        studyFolder,
                        fcsCol="ZBXFN",
                        assay=c("FCM", "CyTOF")){
  cat("Parsing the meta-data from ImmPort\n")
  metaData=unique(metaData[,fcsCol,drop=FALSE])

  if(assay=="FCM"){
    fcs_files=paste(studyFolder,"ResultFiles/Flow_cytometry_result",metaData[,fcsCol],sep="/")
  }else{
    fcs_files=paste(studyFolder,"ResultFiles/CyTOF_result",metaData[,fcsCol],sep="/")
  }
  fsc_markers=sapply(fcs_files,function(x){
    tryCatch({
      fcs=flowCore::read.FCS(x,truncate_max_range=FALSE)
      panel=paste(flowCore::pData(flowCore::parameters(fcs))$name,
                  flowCore::pData(flowCore::parameters(fcs))$desc, sep="-")
      return(paste(panel,collapse="|"))
    },error=function(e){return("DoesNotExist")})
  })
  study=sapply(fsc_markers,function(x){
    if(x=="DoesNotExist" ){return("DoesNotExist")}
    w=which(unique(fsc_markers)==x);
    w=paste0(gsub("/","_",studyFolder),"_",assay,"-",w)
    return(w)
  })
  names(study)=NULL
  inputMeta=data.frame("fcs_files"=fcs_files,"study_id"=study)

  inputMeta=subset(inputMeta,!grepl("DoesNotExist",inputMeta$study_id))
  inputMeta=unique(inputMeta)
  inputMeta$fcs_files=as.character(inputMeta$fcs_files)
  return(inputMeta)
}
