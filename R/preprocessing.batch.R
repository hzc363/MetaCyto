#' Preprocessing fcs files from different studies in batch.
#'
#' It transform and compensate for the raw fcs files and write out the processed
#' data to a new set of fcs files.
#' @param inputMeta A data frame containing 2 columns: a column called
#'   "fcs_files" that contains the location (relative to the working directory)
#'   of each fcs file on the hard drive and a column called "study_id" that
#'   specify what study each fcs file belongs to.
#' @param assay A vector, each element is either "FCM" or "CyTOF" to indicate
#'   the type of cytometry data.
#' @param outpath A string indicating the directory the pre-processed fcs files
#'   will be written to.
#' @param b A positive number used to specify the arcsinh transformation. f(x) =
#'   asinh (b*x) where x is the original value and f(x) is the value after
#'   transformation. The suggested value is 1/150 for flow cytometry (FCM) data
#'   and 1/8 for CyTOF data.
#' @param fileSampleSize An integer specifying the number of events sampled from
#'   each fcs file. If NULL, all the events will be pre-processed and wrote out
#'   to the new fcs files.
#' @param excludeTransformParameters A vector specifying the name of parameters
#'   not to be transformed (left at linear scale).
#' @return Does not return anything. The output is written to the directory
#'   specified by the "outpath".
#' @details The function takes a data frame which specifies the location of the
#'   fcs files and the panels the fcs files belong to. It transform the
#'   cytometry data using the arcsinh transformation. For flow cytometry data, it
#'   compensate the data using the compensation matrix supplied in the fcs file.
#'   the preprocessed fcs files and a table called
#'   "processed_sample_summary.csv" will be wrote out to outpath as well.
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

#' @export

preprocessing.batch = function(inputMeta,
                               assay=c("FCM", "CyTOF"),
                               outpath,
                               b=1/150,
                               fileSampleSize=5000,
                               excludeTransformParameters=c("FSC-A","FSC-W","FSC-H","Time","Cell_length")){

  fcs_param=NULL
  fcs_names= gsub(".*/","",inputMeta$fcs_files)
  excludeTransformParameters=paste(excludeTransformParameters,collapse="|")
  dir.create(outpath,recursive=TRUE)
  if(length(b)==1){b=rep(b,nrow(inputMeta))}
  if(length(assay)==1){assay=rep(assay,nrow(inputMeta))}
  for(std in unique(inputMeta$study_id)){

    cat("Study ID = ",std, " ")

    # 1) identify sample :
    fcs_files=subset(inputMeta$fcs_files,inputMeta$study_id==std)
    fcs=preprocessing(fcsFiles=fcs_files,
                      assay=subset(assay,inputMeta$study_id==std)[1],
                      b=subset(b,inputMeta$study_id==std)[1],
                      fileSampleSize=fileSampleSize,
                      excludeTransformParameters=excludeTransformParameters)

    # 5) outputting resuts
    flowCore::write.FCS(fcs,filename=paste0(outpath,'/',std,".fcs"))

    antibodies=markerFinder(fcs)
    antibodies=paste(antibodies,collapse="|")
    fcs_param=rbind(fcs_param,data.frame("fcs_files"=fcs_files,"study_id"=std,"antibodies"=antibodies))
  }
  fcs_param=cbind("fcs_names"=fcs_names,fcs_param)
  write.csv(fcs_param,paste0(outpath,"/processed_sample_summary.csv"),row.names=FALSE)
  cat("Preprocess result stored in the folder:",outpath, "\n")
}#end of function
