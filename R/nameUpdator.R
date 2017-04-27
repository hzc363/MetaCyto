#' Used to update marker names
#'
#' A function that updates marker names in the files output by the preprocessing.batch function.
#' @param oldNames a vector of marker names you wish to change
#' @param newNames a vector of marker names you wish each oldNames to be changed to.
#' @param files a list of "processed_sample_summary.csv" files output by the preprocessing.batch function, in which the name change will occur.
#' @examples
#' #get meta-data
#' fn=system.file("extdata","fcs_info.csv",package="MetaCyto")
#' fcs_info=read.csv(fn,stringsAsFactors=FALSE,check.names=FALSE)
#' fcs_info$fcs_files=system.file("extdata",fcs_info$fcs_files,package="MetaCyto")
#' # make sure the transformation parameter "b" and the "assay" argument are correct of FCM and CyTOF files
#' b=assay=rep(NA,nrow(fcs_info))
#' b[grepl("CyTOF",fcs_info$study_id)]=1/8
#' b[grepl("FCM",fcs_info$study_id)]=1/150
#' assay[grepl("CyTOF",fcs_info$study_id)]="CyTOF"
#' assay[grepl("FCM",fcs_info$study_id)]="FCM"
#' # preprocessing
#' preprocessing.batch(inputMeta=fcs_info,
#'                     assay=assay,
#'                     b=b,
#'                     outpath="Example_Result/preprocess_output",
#'                     excludeTransformParameters=c("FSC-A","FSC-W","FSC-H","Time","Cell_length"))
#' # Make sure marker names are consistant in different studies
#' files=list.files("Example_Result",pattern="processed_sample",recursive=TRUE,full.names=TRUE)
#' nameUpdator("CD8B","CD8",files)
#' @return Null
#' @export
nameUpdator=function(oldNames,newNames,files){
  oldNames=as.vector(oldNames)
  newNames=as.vector(newNames)
  for(fn in files){
    TB=read.csv(fn,stringsAsFactors=FALSE,check.names=FALSE)
    for(i in 1:length(newNames)){
      TB$antibodies=sapply(TB$antibodies,function(x){
        gsub(oldNames[i],newNames[i],x,fixed=TRUE)
      })
    }
    write.csv(TB,fn,row.names=FALSE)
  }
}
