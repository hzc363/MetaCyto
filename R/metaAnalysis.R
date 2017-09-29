#' Performs meta-analysis
#'
#' A function that performs meta-analysis
#' @param value A string to specify the column name of the dependent variable
#'   (y)
#' @param variableOfInterst A string to specify the column name of the
#'   independent variable of interest (x1)
#' @param otherVariables A string vector to specify the column names of
#'   independent variables included in the regression model other than the
#'   variableOfInterst.
#' @param studyID A string to specify the column name of study ID.
#' @param data A data frame containing the data
#' @param CILevel A number between 0 to 1, used to specify the confidence
#'   interval to be plotted in the forest plot.
#' @param main A string to specify the title of the forest plot
#' @param ifScale A vector of two logic values, specifying if the dependent
#'   variable and the variableOfInterst should be scaled when calculating the
#'   effect size.
#' @param cex A number specifying the amount by which plotting text and symbols
#'   should be scaled relative to the default in the forest plot.
#' @return Returns data frame describing the effect size of variableOfInterst on
#'   value in each individual studies, as well as the over all effect size.
#' @examples
#' library(dplyr)
#' #collect all summary statistics
#' fn=system.file("extdata","",package="MetaCyto")
#' files=list.files(fn,pattern="cluster_stats_in_each_sample",recursive=TRUE,
#'                  full.names=TRUE)
#' fcs_stats=collectData(files,longform=TRUE)
#' # Collect sample information
#' files=list.files(fn,pattern="sample_info",recursive=TRUE,full.names=TRUE)
#' sample_info=collectData(files,longform=FALSE)
#' # join the cluster summary statistics with sample information
#' all_data=inner_join(fcs_stats,sample_info,by="fcs_files")
#'
#' # plot forrest plot to see if the proportion of CCR7+ CD8 T cell
#' # is affected by age (while controlling for GENDER)
#' L="CD3+|CD4-|CD8+|CCR7+"
#' dat=subset(all_data,all_data$parameter_name=="fraction"&
#'             all_data$label==L)
#' MA=metaAnalysis(value="value",variableOfInterst="SUBJECT_AGE",main=L,
#'                 otherVariables=c("GENDER"),studyID="study_id",
#'                 data=dat,CILevel=0.95,ifScale=c(TRUE,FALSE))
#' @importFrom metafor rma.uni forest
#' @export
metaAnalysis=function(value,variableOfInterst,otherVariables,
                      studyID,data,CILevel,main,
                      ifScale=c(TRUE,FALSE),cex=1){
  study_result=NULL
  for(std in unique(data[,studyID])){
    sub_data=subset(data,data[,studyID]==std)
    sub_data = sub_data[,c(value,variableOfInterst,otherVariables)]
    sub_data=na.omit(sub_data)
    NL=apply(sub_data[,c(value,variableOfInterst,otherVariables)],2,function(x){length(unique(x))})
    if(!all(NL>1)){cat(std,"is skipped. One of the variable have 0 variance.\n");next}
    if(ifScale[1]){sub_data[,value]=scale(sub_data[,value])}
    if(ifScale[2]){sub_data[,variableOfInterst]=scale(sub_data[,variableOfInterst])}
    colnames(sub_data)[1]="Y"
    LM=lm(Y~.,data=sub_data)
    CE=t(summary(LM)$coefficients[2,])
    CI=confint(LM,parm=2,level=CILevel)
    t1=cbind("study_id"=std,data.frame(CE,check.names=FALSE),
             data.frame(CI,check.names=FALSE),"N"=nrow(sub_data))
    study_result=rbind(study_result,t1)
  }
  study_result=na.omit(study_result)
  res=metafor::rma.uni(yi=study_result$Estimate, vi=(study_result$`Std. Error`)^2)
  metafor::forest(res, slab=study_result$study_id,main=main,
                  xlab="Effect Size", mlab="RE Model for All Studies",cex=cex)
  t1=data.frame("Summary",res$b[1],res$se,res$zval,res$pval,res$ci.lb,res$ci.ub,
                sum(study_result$N))
  names(t1)=names(study_result)
  study_result=rbind(study_result,t1)
  rownames(study_result)=NULL
  return(study_result)
}
