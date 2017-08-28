#' Derive summary statistics for clusters
#'
#' A function that derives summary statistics for clusters.
#' @param fcsFrame A flow Frame object returned from preprocessing function.
#'   Must contain a parameter called "sample_id" that specify the origin of
#'   each cell using integer IDs.
#' @param clusterList A list, each element should be a vector containing the IDs
#'   of all cells that belong to a cluster.
#' @param fcsNames A vector of fcs file names. Each element corresponds to an
#'   integer ID in the "sample_id" parameter in fcsFrame.
#' @return Returns a data frame, rows correspond to each fcs file, columns
#'   correspond to MFI of markers or fractions.
#' @examples
#' # Find fcs files
#' files=system.file("extdata","SDY420/ResultFiles/CyTOF_result",
#'                   package="MetaCyto")
#' files=list.files(files,pattern="fcs$",full.names=TRUE)
#' # Preprocess
#' fcs = preprocessing(fcsFiles=files,assay ="CyTOF",b=1/8)
#' # Search clusters
#' cluster_list=searchCluster(fcsFrame=fcs,
#'                            clusterLabel=c("CD3+|CD8+","CD3-|CD19+"))
#' # derive summary statistics
#' cluster_stats=clusterStats(fcsFrame=fcs,
#'                            clusterList=cluster_list$clusterList,
#'                            fcsNames=files)
#' @importFrom flowCore exprs
#' @export
clusterStats=function(fcsFrame,clusterList,fcsNames){
  antibodies=markerFinder(fcsFrame)
  expr=flowCore::exprs(fcsFrame)
  colnames(expr)=antibodies
  CL=clusterList
  CL_label=names(clusterList)
  sample_name=fcsNames
  expr_sample=expr[,"SAMPLE_ID"]
  expr=expr[,colnames(expr)!="SAMPLE_ID"]
  antibodies=antibodies[antibodies!="SAMPLE_ID"]
  SP_stat=NULL
  for(i in 1:length(CL)){
    for(j in unique(expr_sample)){
      fileSampleSize=sum(expr_sample==j)
      w=intersect( which(expr_sample==j), CL[[i]] )
      frac=length(w)/fileSampleSize
      if(length(w)>1){M=apply(expr[w,],2,median)}
      if(length(w)==1){M=expr[w,]}
      if(length(w)==0){M=rep(NA,ncol(expr))}
      t1=c(paste0("cluster",i),CL_label[i],sample_name[j],M,frac)
      SP_stat=rbind(SP_stat,t1)
    }
  }
  colnames(SP_stat)=c("cluster_id","label","fcs_names",antibodies,"fraction")
  rownames(SP_stat)=NULL
  SP_stat=data.frame(SP_stat,stringsAsFactors = FALSE,check.names = FALSE)
  SP_stat[,1:3]=lapply(SP_stat[,1:3],as.character)
  SP_stat[,-c(1:3)]=lapply(SP_stat[,-c(1:3)],as.numeric)

  return(SP_stat)
}
