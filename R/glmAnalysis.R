#' Perform generalized linear model analysis to estimate effect size.
#'
#' A function that performs generalized linear model analysis to estimate effect size.
#' @param value A string to specify the column name of the dependent variable (y)
#' @param variableOfInterst A string to specify the column name of the independent variable of interest (x1)
#' @param otherVariables A string vector to specify the column names of independent variables included in the regression model other than the variableOfInterst.
#' @param parameter A string to specify what summary statistics is the dependent variable.
#' @param studyID A string to specify the column name of study ID.
#' @param label A string to specify the name the column that contains the cluster label or name.
#' @param data A data frame containing the data. Usually a long form data frame returned by collectData.
#' @param CILevel A number between 0 to 1, used to specify the confidence interval to be plotted in the forest plot.
#' @param ifScale A vector of two logic values, specifying if the dependent variable and the variableOfInterst should be scaled when calculating the effect size.
#' @return Returns data frame describing the overall effect size of variableOfInterst on value.  May be slightly different from the value reported from the function metaAnalysis.
#' @details The function use the model value ~ variableOfInterst + otherVariables + studyID to estimate the effect size. Use it as a screening tool. Use metaAnalysis function to analyze an effect size in more detail.
#' @examples
#' #collect all summary statistics
#'fn=system.file("extdata","",package="MetaCyto")
#'files=list.files(fn,pattern="cluster_stats_in_each_sample",recursive=TRUE,full.names=TRUE)
#'fcs_stats=collectData(files,longform=TRUE)
#'# Collect sample information
#'files=list.files(fn,pattern="sample_info",recursive=TRUE,full.names=TRUE)
#'sample_info=collectData(files,longform=FALSE)
#'# join the cluster summary statistics with sample information
#'all_data=inner_join(fcs_stats,sample_info,by="fcs_files")
#'# See the fraction of what clusters are affected by age (while controlling for GENDER)
#'GA=glmAnalysis(value="value",variableOfInterst="SUBJECT_AGE",parameter="fraction",
#'               otherVariables=c("GENDER"),studyID="study_id",label="label",
#'               data=all_data,CILevel=0.95,ifScale=c(TRUE,FALSE))
#' @export
glmAnalysis=function(value="value",variableOfInterst="SUBJECT_AGE",parameter,
                     otherVariables=c("GENDER"),studyID="study",label="label",
                     data,CILevel=0.95,ifScale=c(TRUE,FALSE)){
  result=NULL
  for(L in unique(data[,label])){
    sub_data=data[data$parameter_name==parameter&data[,label]==L,]
    sub_data=na.omit(sub_data)
    if(length(unique(sub_data[,studyID]))>1){
      sub_data=sub_data[,c(value,variableOfInterst,otherVariables,studyID)]
    }else{sub_data=sub_data[,c(value,variableOfInterst,otherVariables)]}

    NL=apply(sub_data,2,function(x){length(unique(x))})
    if(!all(NL>1)){cat(L,"is skipped. One of the variable have 0 variance.\n");next}
    if(ifScale[1]){sub_data[,value]=scale(sub_data[,value])}
    if(ifScale[2]){sub_data[,variableOfInterst]=scale(sub_data[,variableOfInterst])}
    colnames(sub_data)[1]="Y"
    LM=lm(Y~.,data=sub_data)
    CE=t(summary(LM)$coefficients[2,])
    CI=confint(LM,parm=2,level=CILevel)
    t1=cbind("label"=L,data.frame(CE,check.names=FALSE),
             data.frame(CI,check.names=FALSE),"N"=nrow(sub_data))
    result=rbind(result,t1)
  }
  if(!is.null(result)){
    rownames(result)=NULL
    colnames(result)=c("label","Effect_size","SE","t_value","p_value","lower","upper","N")
  }

  return(result)
}
