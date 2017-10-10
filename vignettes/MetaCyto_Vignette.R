## ---- echo=FALSE, results='asis'-----------------------------------------
fn=system.file("extdata","fcs_info.csv",package="MetaCyto")
fcs_info=read.csv(fn,stringsAsFactors=FALSE,check.names=FALSE)
fcs_info=fcs_info[order(fcs_info$fcs_files),]
row.names(fcs_info)=NULL
knitr::kable(head(fcs_info, 6))

## ---- echo=FALSE, results='asis'-----------------------------------------
fn=system.file("extdata","sample_info_vignette.csv",package="MetaCyto")
sample_info=read.csv(fn,stringsAsFactors=FALSE,check.names=FALSE)
sample_info=sample_info[order(sample_info$fcs_files),]
row.names(sample_info)=NULL
knitr::kable(head(sample_info, 6))

## ---- fig.show='hold'----------------------------------------------------
# allocate variables
b=assay=rep(NA,nrow(fcs_info))

# define transformation parameter for each fcs files
b[grepl("CyTOF",fcs_info$study_id)]=1/8
b[grepl("FCM",fcs_info$study_id)]=1/150

# define cytometry type for each fcs files
assay[grepl("CyTOF",fcs_info$study_id)]="CyTOF"
assay[grepl("FCM",fcs_info$study_id)]="FCM"

## ---- results="hide", message=FALSE, warning=FALSE-----------------------
library(MetaCyto)

# find example data in the MetaCyto package. You won't need this line when running your actual meta-analysis
fcs_info$fcs_files=system.file("extdata",fcs_info$fcs_files,package="MetaCyto") 

# preprocessing the data
preprocessing.batch(inputMeta=fcs_info,
                     assay=assay,
                     b=b,
                     outpath="Example_Result/preprocess_output",
                     excludeTransformParameters=c("FSC-A","FSC-W","FSC-H","Time","Cell_length"))

## ----results='asis',echo=TRUE--------------------------------------------
# collect preprocessing information
files=list.files("Example_Result",pattern="processed_sample",recursive=TRUE,full.names=TRUE)
panel_info=collectData(files,longform=FALSE)

# analyze the panels
PS=panelSummary(panelInfo = panel_info,
                folder = "Example_Result",
                cluster=FALSE,
                width=30,
                height=20)
knitr::kable(head(PS))

## ----results='asis',echo=TRUE--------------------------------------------
sort(rownames(PS))

## ----results='asis',echo=TRUE--------------------------------------------
nameUpdator(oldNames=c("CD8B"), newNames=c("CD8"), files=files)

## ---- results="hide", message=FALSE, warning=FALSE-----------------------
#define parameters that we don't want to cluster
excludeClusterParameters=c("FSC-A","FSC-W","FSC-H","SSC-A","SSC-W","SSC-H","Time",
                           "CELL_LENGTH","DEAD","DNA1","DNA2")

# Find and label clusters in the data. The default cluster functions 
# are FlowSOM.MC. Here, flowHC, a hiarachical clustering function is used instead
# to show how non-default functions can be used. 
cluster_label=autoCluster.batch(preprocessOutputFolder="Example_Result/preprocess_output",
                                 excludeClusterParameters=excludeClusterParameters,
                                 labelQuantile=0.95,
                                 clusterFunction=flowHC)

## ---- results="hide", message=FALSE, warning=FALSE-----------------------
cluster_label=c(cluster_label,"CD3+|CD4-|CD8+|CCR7+")

## ---- results="hide", message=FALSE, warning=FALSE-----------------------
searchCluster.batch(preprocessOutputFolder="Example_Result/preprocess_output",
              outpath="Example_Result/search_output",
              clusterLabel=cluster_label)

## ---- out.width = "700px",echo=FALSE-------------------------------------
fn=system.file("extdata","density_plot.png",package="MetaCyto")
knitr::include_graphics(fn)

## ---- results="hide", message=FALSE, warning=FALSE-----------------------
library(dplyr)
# Collect Summary statistics generated in step 3
files=list.files("Example_Result/search_output",pattern="cluster_stats_in_each_sample",recursive=TRUE,full.names=TRUE)
fcs_stats=collectData(files,longform=TRUE)

# Get sample information generated in step 1
fn=system.file("extdata","sample_info_vignette.csv",package="MetaCyto")
sample_info=read.csv(fn,stringsAsFactors=FALSE,check.names=FALSE)

# find data in the MetaCyto package. You won't need this line when running your actual meta-analysis
sample_info$fcs_files=system.file("extdata",sample_info$fcs_files,package="MetaCyto") 

# join the cluster summary statistics with sample information
all_data=inner_join(fcs_stats,sample_info,by="fcs_files")

## ---- fig.show='hold', message=FALSE, warning=FALSE , fig.height=4, fig.width=6----
# See the fraction of what clusters are affected by age (while controlling for GENDER)
GA=glmAnalysis(value="value",variableOfInterst="SUBJECT_AGE",parameter="fraction",
               otherVariables=c("GENDER"),studyID="study_id",label="label",
               data=all_data,CILevel=0.95,ifScale=c(TRUE,FALSE))
GA=GA[order(GA$Effect_size),]

# To save space, only cell populations with short cell definitions are plotted
GA$label=as.character(GA$label)
w = which(nchar(GA$label)<30)
GA = GA[w,]

# plot the results
plotGA(GA)

## ---- fig.show='hold', message=FALSE, warning=FALSE , fig.height=4, fig.width=6----
# Subset the data to only include effect size of age on the proportion of "CD3+ CD4- CD8+ CCR7+"
L="CD3+|CD4-|CD8+|CCR7+"
dat=subset(all_data,all_data$parameter_name=="fraction"&
             all_data$label==L)

# run the metaAnalysis function
MA=metaAnalysis(value="value",variableOfInterst="SUBJECT_AGE",main=L,
                  otherVariables=c("GENDER"),studyID="study_id",
                  data=dat,CILevel=0.95,ifScale=c(TRUE,FALSE))

## ---- echo=FALSE---------------------------------------------------------
unlink("Example_Result",recursive = TRUE)

## ---- fig.show='hold', results="hide", message=FALSE, warning=FALSE , fig.height=4, fig.width=6----

# read meta-data of SDY736
fn=system.file("extdata","SDY736/SDY736-DR19_Subject_2_Flow_cytometry_result.txt",package="MetaCyto")
meta_data=read.table(fn,sep='\t',header=TRUE)

# Organize fcs file into panels
fn=system.file("extdata","SDY736",package="MetaCyto")
fcs_info_SDY736=fcsInfoParser(metaData=meta_data,studyFolder=fn,
                               fcsCol="FILE_NAME",assay="FCM")

## ---- fig.show='hold', results="hide", message=FALSE, warning=FALSE , fig.height=4, fig.width=6----
# read meta-data of SDY420
fn=system.file("extdata","SDY420/SDY420-DR19_Subject_2_CyTOF_result.txt",package="MetaCyto")
meta_data=read.table(fn,sep='\t',header=TRUE)

# Organize fcs file into panels
fn=system.file("extdata","SDY420",package="MetaCyto")
fcs_info_SDY420=fcsInfoParser(metaData=meta_data,studyFolder=fn,
                               fcsCol="FILE_NAME",assay="CyTOF")


## ---- fig.show='hold',  results="hide",message=FALSE, warning=FALSE , fig.height=4, fig.width=6----
# Combine fcs info
fcs_info=rbind(fcs_info_SDY420,fcs_info_SDY736)

## ---- fig.show='hold', results="hide", message=FALSE, warning=FALSE , fig.height=4, fig.width=6----
# read meta-data of SDY736
fn=system.file("extdata","SDY736/SDY736-DR19_Subject_2_Flow_cytometry_result.txt",package="MetaCyto")
meta_data=read.table(fn,sep='\t',header=TRUE)

# Find the AGE, GENDER info from selected_data
fn=system.file("extdata","SDY736",package="MetaCyto")
sample_info_SDY736=sampleInfoParser(metaData=meta_data,
                                   studyFolder=fn,
                                   assay="FCM",
                                   fcsCol="FILE_NAME",
                                   attrCol=c("SUBJECT_AGE","GENDER"))

## ---- fig.show='hold', results="hide", message=FALSE, warning=FALSE , fig.height=4, fig.width=6----
# read meta-data of SDY420
fn=system.file("extdata","SDY420/SDY420-DR19_Subject_2_CyTOF_result.txt",package="MetaCyto")
meta_data=read.table(fn,sep='\t',header=TRUE)

# Find the AGE, GENDER info from selected_data
fn=system.file("extdata","SDY736",package="MetaCyto")
sample_info_SDY420=sampleInfoParser(metaData=meta_data,
                                   studyFolder=fn,
                                   assay="CyTOF",
                                   fcsCol="FILE_NAME",
                                   attrCol=c("SUBJECT_AGE","GENDER"))


## ---- fig.show='hold',  results="hide",message=FALSE, warning=FALSE , fig.height=4, fig.width=6----
# Combine sample info
sample_info=rbind(sample_info_SDY420,sample_info_SDY736)

