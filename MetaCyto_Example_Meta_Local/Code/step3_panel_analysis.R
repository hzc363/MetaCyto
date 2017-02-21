library(MetaCyto)
library(dplyr)

# The purpose of this step is to make sure that the marker names are consistent between 
# studies. This is immportant for combining data later. 

#########################################
########## Examine the Panels ##########
########################################

# collect data output from the step1. 
files=list.files("Result",pattern="processed_sample",recursive=T,full.names=T)
panel_info=collectData(files,longform=F)

# summarize markers in all markers. The panel summary is written to "Result" Folder
PS=panelSummary(panel_info,"Result",cluster=F,width=30,height=20) 


###################################
########## Updata Names ##########
##################################

# CD8 is called "CD8B" in one study and is called "CD8" in another. 
# Here, we change CD8B to CD8 for consistency.
ab_names=sort(rownames(PS))
nameUpdator("CD8B","CD8",files)


#####################################################
########## Examine the Panels after update ##########
#####################################################
panel_info=collectData(files,longform=F)
PS=panelSummary(panel_info,"Result",cluster=F,width=30,height=20)



