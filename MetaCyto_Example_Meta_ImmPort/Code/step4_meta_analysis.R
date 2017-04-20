library(MetaCyto)
library(dplyr)

###############################################
########## Collect Summary statistics  ########
###############################################
# Collect Summary statistics generated in step 3
files=list.files("Result/search_output",pattern="cluster_stats_in_each_sample",recursive=T,full.names=T)
fcs_stats=collectData(files,longform=T)


#######################################################
########## Collect collect sample information  ########
#######################################################

# Collect sample information generated in step 1
files=list.files("Result",pattern="sample_info",recursive=T,full.names=T)
sample_info=collectData(files,longform=F)

# join the cluster summary statistics with sample information
all_data=inner_join(fcs_stats,sample_info,by="fcs_files")



####################################################################
########## Estimation effect size  of age for all clusters  ########
####################################################################

# See the fraction of what clusters are affected by age (while controlling for GENDER)
GA=glmAnalysis(value="value",variableOfInterst="SUBJECT_AGE",parameter="fraction",
               otherVariables=c("GENDER"),studyID="study_id",label="label",
               data=all_data,CILevel=0.95,ifScale=c(T,F))
GA=GA[order(GA$Effect_size),]

# plot the results
plotGA(GA)


##############################################################
########## Analyze on cluster  in detail  ####################
##############################################################
# Let's take a closer look to see if CCR7+ CD8 T cell is affected by age (while controlling for GENDER)
L="CD3+|CD4-|CD8+|CCR7+"
dat=subset(all_data,all_data$parameter_name=="fraction"&
             all_data$label==L)
MA=metaAnalysis(value="value",variableOfInterst="SUBJECT_AGE",main=L,
                  otherVariables=c("GENDER"),studyID="study_id",
                  data=dat,CILevel=0.95,ifScale=c(T,F))

