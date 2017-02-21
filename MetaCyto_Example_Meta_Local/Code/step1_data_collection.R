library(dplyr)
library(tidyr)
library(MetaCyto)

# This example is a light weighted example of how to use MetaCyto to joint analyze your own local data. 

###################################################
########## Check working directory first ##########
###################################################
wd=getwd()
if(grepl("MetaCyto_Example_Meta_Local$",wd)==F){print("Please change the working directory to 'MetaCyto_Example'")}


#############################################
########## Collect data Manually  ##########
############################################
# For your own data, the data collection step need to be done manually.
# First, you need to prepare a fcs_info table like the one below. 
#  The fcs_files column specifies the location of the fcs files.
# The study_id column specifies that panels. 
#Each panel is a collection of fcs files that have the same markers and are from the same study. 

fcs_info=read.csv("Data/fcs_info.csv")
print(fcs_info)


# If you want to associate fcs data with some phenotype or genotype data, such as age, gender or SNP, 
# you need to collect the sample information as well, like the one below. 
#  The fcs_files column specifies the location of the fcs files.
# Other two columns specifies the age and gender of each subject corresponding to each fcs file. 
# you can include whatever sample information in the table. 
sample_info=read.csv("Data/sample_info.csv")
print(sample_info)


