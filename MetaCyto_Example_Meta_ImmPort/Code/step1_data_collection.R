library(dplyr)
library(tidyr)
library(MetaCyto)

# This example is a light weighted example of how to use MetaCyto to analyze data downloaded from 
# ImmPort. First, we collect fcs files from SDY736(flow cytometry data) and SDY420(CyTOF data). 

###################################################
########## Check working directory first ##########
###################################################
wd=getwd()
if(grepl("MetaCyto_Example_Meta_ImmPort$",wd)==F){print("Please change the working directory to 'MetaCyto_Example'")}


################################################
########## SDY736 FCM data collection ##########
################################################

# Read Meta-data from ImmPort
meta_data=read.table("Data/SDY736/SDY736-DR19_Subject_2_Flow_cytometry_result.txt",sep='\t',header=T)

# Organize fcs file into panels
fcs_info_SDY736=fcsInfoParser(metaData=meta_data,
                       studyFolder="Data/SDY736",
                       fcsCol="FILE_NAME",
                       assay="FCM")

# Find the AGE, GENDER info from selected_data
sample_info_SDY736=sampleInfoParser(metaData=meta_data,
                             studyFolder="Data/SDY736",
                             assay="FCM",
                             fcsCol="FILE_NAME",
                             attrCol=c("SUBJECT_AGE","GENDER"))


##################################################
########## SDY420 CyTOF data collection ##########
##################################################

meta_data=read.table("Data/SDY420/SDY420-DR19_Subject_2_CyTOF_result.txt",sep='\t',header=T)

fcs_info_SDY420=fcsInfoParser(metaData=meta_data,
                       studyFolder="Data/SDY420",
                       fcsCol="FILE_NAME",
                       assay="CyTOF")

sample_info_SDY420=sampleInfoParser(metaData=meta_data,
                             studyFolder="Data/SDY420",
                             assay="CyTOF",
                             fcsCol="FILE_NAME",
                             attrCol=c("SUBJECT_AGE","GENDER"))

#########################################################
########## Combine Data from SDY420 and SDY736 ##########
#########################################################
sample_info=rbind(sample_info_SDY736,sample_info_SDY420)
fcs_info=rbind(fcs_info_SDY736,fcs_info_SDY420)
dir.create("Result")
write.csv(sample_info,"Result/sample_info.csv",row.names=F)
write.csv(fcs_info,"Result/fcs_info.csv",row.names=F)

#################################
########## References ##########
################################
# Data used in the Example is a subset of data from SDY420[1] and SDY736[2] on ImmPort[3]. 
# 
# 
# [1] Whiting, Chan C., et al. "Large-scale and comprehensive immune profiling and functional analysis of normal human aging." PloS one 10.7 (2015): e0133627.
# 
# [2] Wertheimer, Anne M., et al. "Aging and cytomegalovirus infection differentially and jointly affect distinct circulating T cell subsets in humans." The Journal of Immunology 192.5 (2014): 2143-2155.
# 
# [3]Bhattacharya, Sanchita, et al. "ImmPort: disseminating data to the public for the future of immunology." Immunologic research 58.2-3 (2014): 234-239.

