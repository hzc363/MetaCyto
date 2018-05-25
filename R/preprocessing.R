#' Preprocess fcs files from a single experiment
#'
#' A function that preprocesses fcs files from a single experiment.
#' @param fcsFiles A vector specifying the location of all fcs files.
#' @param assay Either "FCM" or "CyTOF" to indicate the type of cytometry data.
#' @param b A positive number used to specify the arcsinh transformation. f(x) =
#'   asinh (b*x) where x is the original value and f(x) is the value after
#'   transformation. The suggested value is 1/150 for flow cytometry (FCM) data
#'   and 1/8 for CyTOF data. If b = 0, the transformation is skipped.
#' @param fileSampleSize An integer specifying the number of events sampled from
#'   each fcs file. If NULL, all the events will be pre-processed and wrote out
#'   to the new fcs files.
#' @param excludeTransformParameters A vector specifying the name of parameters
#'   not to be transformed (left at linear scale).
#' @param compFiles A vector specifying the paths of user supplied compensation
#'   matrix for each fcs file. The matrix must be stored in csv files.
#' @return Returns a flowFrame object containing the preprocessed cytometry
#'   data. Cells from different fcs files are combined into one flow frame. A
#'   new parameter, sample_id, is introduced to indicate the origin of each
#'   cell.
#' @examples
#' # Find fcs files
#' files=system.file("extdata","SDY420/ResultFiles/CyTOF_result",
#'                   package="MetaCyto")
#' files=list.files(files,pattern="fcs$",full.names=TRUE)
#' # Preprocess
#' fcs = preprocessing(fcsFiles=files,assay ="CyTOF",b=1/8)
#' @importFrom flowCore read.flowSet fsApply keyword logTransform transformList
#'   transform compensate arcsinhTransform colnames
#' @export
preprocessing=function(fcsFiles,
                       assay=c("FCM", "CyTOF"),
                       b=1/200,
                       fileSampleSize=5000,
                       compFiles=NULL,
                       excludeTransformParameters=c("FSC-A","FSC-W","FSC-H","Time","Cell_length")){

  fcs_param=NULL
  fcs_names= gsub(".*/","",fcsFiles)
  excludeTransformParameters=paste(excludeTransformParameters,collapse="|")


  cat("Preprocessing \n")

  # 1) identify sample :
  fcsFiles=as.character(fcsFiles)

  # 2) read the fcs files
  fcs = flowCore::read.flowSet(fcsFiles,transformation="linearize",alter.names=FALSE,truncate_max_range=FALSE)
  if(!is.null(fileSampleSize)& (!is.na(fileSampleSize)) ){
    fcs = flowCore::fsApply(fcs, function(f){
      L=nrow(f@exprs)
      if(L>fileSampleSize){
        f@exprs=f@exprs[sample(1:L,fileSampleSize),]
      }
      return(f)
    })
  }

  # 3) compensation and transformation
  if (assay == "FCM") {
    # retreiving fluorescent biomarkers
    w=which(!grepl(excludeTransformParameters,flowCore::colnames(fcs),ignore.case = TRUE))
    biomarker_vector = flowCore::colnames(fcs)[w]
    # identify the fcs file format
    version = flowCore::fsApply(fcs, function(frame) {
      return(flowCore::keyword(frame, "FCSversion")[[1]])
    })
    unique_version = as.numeric(unique(as.vector(version)))

    # compensation, transformation
    ## if the file version is 2.0, then code will log transform the data
    if (unique_version == 2) {
      trans = flowCore::logTransform()
      translist = flowCore::transformList(biomarker_vector,trans)
      fcs = flowCore::transform(fcs, translist)
    } else if (unique_version == 3) {
      # check if user have provided the compmatrix
      if(!is.null(compFiles)){
        compList=lapply(compFiles,function(x){
          if(is.na(x)){return(as.matrix(NA))}else{return(as.matrix(read.csv(x,row.names = 1,check.names = FALSE)))}
          })
        for(id in 1:length(fcs)){
          if(!is.na(compList[[id]][1,1])){
            fcs[[id]]@description$SPILL <- compList[[id]]
            }
        }
      }

      ### function first checks if the spill matrix in the flowframe is an actual spill matrix or an identity matrix
      if (is.null(flowCore::keyword(fcs[[1]], "SPILL")[[1]]) == FALSE) {
        check = flowCore::fsApply(fcs, function(x) {
          result = isSymmetric(flowCore::keyword(x, "SPILL")[[1]])
        })
        ### if all flowframes have spill matrices then my function
        ### will apply them to each flowframe in the flowset for compensation
        if (unique(check)[[1]] == FALSE) {
          fcs = flowCore::fsApply(fcs, function(x) {
            new_frame = flowCore::compensate(x, flowCore::keyword(x, "SPILL")[[1]])
            return(new_frame)
          })
        }
      }
      if(b != 0){
        trans = flowCore::arcsinhTransform(transformationId="defaultArcsinhTransform",a=0,b=b,c=0)
        translist = flowCore::transformList(biomarker_vector, trans)
        fcs = flowCore::transform(fcs, translist)
      }

    }
  } else if (assay == "CyTOF") {
    w=which(!grepl(excludeTransformParameters,flowCore::colnames(fcs),ignore.case = TRUE))
    biomarker_vector = flowCore::colnames(fcs)[w]

    if(b != 0){
    trans = flowCore::arcsinhTransform(transformationId = "defaultArcsinhTransform", a = 0, b = b, c = 0)
    translist = flowCore::transformList(biomarker_vector, trans)
    fcs = flowCore::transform(fcs, translist)
    }
  }
  fcs=set2Frame(fcs)
  return(fcs)
}
