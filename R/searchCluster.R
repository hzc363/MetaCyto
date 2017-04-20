#' search for clusters using pre-defined labels
#'
#' A function that search for clusters using pre-defined labels (cell definitions).
#' @param fcsFrame A flowFrame object.
#' @param clusterLabel a vector of labels, such as "CD3+|CD4+|CD8-". Each marker is followed by "+" or "-" and are sperated by "|".
#' @param cutoff a vector of cutoff values to bisect the distribution of each
#'   marker. The names of the vector should be the same as the marker names. If
#'   NULL, the cutoff value will be determined automatically.
#' @param rmNull True or False. Used to specifiy if a cluster with 0 cells
#'   should be returned or not.
#' @param preGate A character string specifiying the gated used to clean up the data. For example, use "PI-" to only analyze live cell. Or use "Cell_Length+" to only analyze non-debris.
#' @return Returns a list with two components: 1) clusterList is a list in which
#'   each element of the list is a vector containing the ID of all cells in a
#'   cluster. The names corresponds to the labels specified in clusterLabel. 2)
#'   cutoff, contains a vector of cutoff values used to bisect each marker.
#' @examples
#' # Find fcs files
#' files=system.file("extdata","SDY420/ResultFiles/CyTOF_result",package="MetaCyto")
#' files=list.files(files,pattern="fcs$",full.names=TRUE)
#' # Preprocess
#' fcs = preprocessing(fcsFiles=files,assay ="CyTOF",b=1/8)
#' # Search clusters
#' cluster_list=searchCluster(fcsFrame=fcs,
#'                            clusterLabel=c("CD3+|CD8+","CD3-|CD19+"))
#' @export
searchCluster=function(fcsFrame,
                       clusterLabel,
                       cutoff=NULL,
                       rmNull=TRUE,
                       preGate=NULL){

  #prepare exclude parameters
  clusterLabel=toupper(clusterLabel)
  antibodies=markerFinder(fcsFrame)

  # Get expression matrix
  expr=flowCore::exprs(fcsFrame);
  colnames(expr)=antibodies


  # pre-gating
  P1=1:nrow(expr)
  if(!is.null(preGate)){
    preGate=toupper(preGate)
    pre_vec = strsplit(preGate,split="&|\\|")[[1]]
    pre_marker = gsub("\\+$|-$|\\^NE$|\\^LO$|\\^HI$","",pre_vec )

    cutoff0=apply(X=expr[,pre_marker,drop=FALSE],MARGIN=2,FUN=findCutoff)

    P1=sapply(pre_vec,function(x){
      m=gsub("\\+$|-$|\\^NE$|\\^LO$|\\^HI$","",x)
      if(grepl("\\+$",x)){w=which(expr[,m]>cutoff0[m])}
      if(grepl("-$",x)){w=which(expr[,m]<cutoff0[m])}
      if(grepl("\\^NE$",x)){w=which(expr[,m]<triS[2,m])}
      if(grepl("\\^LO$",x)){w=which(expr[,m]<triS[1,m])&expr[,m]>triS[2,m]}
      if(grepl("\\^HI$",x)){w=which(expr[,m]>triS[1,m])}
      return(w)
    })
    if(length(pre_vec)>1){P1=Reduce(intersect,P1)}
  }

  # subset expr
  AB=unlist(strsplit(clusterLabel,split="&|\\|"))
  AB=sapply(AB,function(x){gsub("\\+$|-$|\\^NE$|\\^LO$|\\^HI$","",x)})
  w=which(antibodies%in%AB)
  expr=expr[,w,drop=FALSE]

  #find labels whose markers are included in the study
  CL_label=clusterLabel
  CL_label2=strsplit(clusterLabel,split="&|\\|")
  in_study=sapply(CL_label2,function(x){
    Ab_in_label=gsub("\\+$|-$|\\^NE$|\\^LO$|\\^HI$","",x)
    if(length(setdiff(Ab_in_label,antibodies))>0){return(FALSE)}else{return(TRUE)}
  })
  CL_label=CL_label[in_study]
  CL_label2=CL_label2[in_study]
  if(length(CL_label)<1){return(list())}

  #find cutoff of each parameter
  if(!setequal(AB,names(cutoff))){
    t1=apply(X=expr[P1,],MARGIN=2,FUN=findCutoff)
    if(!is.null(cutoff)){
      cutoff=cutoff[names(cutoff)%in%AB]
      t1[names(cutoff)]=cutoff
      }
    cutoff=t1
  }

  #find each cluster
  CL=list()
  for(i in 1:length(CL_label2)){
    rows=sapply(CL_label2[[i]],function(x){
      m=gsub("\\+$|-$|\\^NE$|\\^LO$|\\^HI$","",x)
      if(grepl("\\+$",x)){w=which(expr[,m]>cutoff[m])}
      if(grepl("-$",x)){w=which(expr[,m]<cutoff[m])}
      if(grepl("\\^NE$",x)){w=which(expr[,m]<triS[2,m])}
      if(grepl("\\^LO$",x)){w=which(expr[,m]<triS[1,m])&expr[,m]>triS[2,m]}
      if(grepl("\\^HI$",x)){w=which(expr[,m]>triS[1,m])}
      return(w)
    })
    if(length(CL_label2[[i]])>1){rows=Reduce(intersect,rows)}
    CL[[i]]=rows
  }
  if(rmNull==TRUE){
    CL_count=sapply(CL,length)
    CL=CL[CL_count>1]
    CL_label=CL_label[CL_count>1]
  }

  if(length(CL)<1){return(list())}
  CL=lapply(CL,intersect,y=P1)

  names(CL)=CL_label
  Result=list("clusterList"=CL,"cutoff"=cutoff)
  return(Result)
}
