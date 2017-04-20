library(MetaCyto)

#################################################
########## Prepare for preprossessing ##########
################################################

# Read the fcs_info file from step 1
fcs_info=read.csv("Result/fcs_info.csv",stringsAsFactors=F,check.names=F)

# make sure the transformation parameter "b" and the "assay" argument are correct of FCM
# and CyTOF files
b=assay=rep(NA,nrow(fcs_info))
b[grepl("CyTOF",fcs_info$study_id)]=1/8
b[grepl("FCM",fcs_info$study_id)]=1/150
assay[grepl("CyTOF",fcs_info$study_id)]="CyTOF"
assay[grepl("FCM",fcs_info$study_id)]="FCM"


#############################################
########## Perform preprossessing ##########
############################################

preprocessing.batch(inputMeta=fcs_info,
                    assay=assay,
                    b=b,
                    outpath="Result/preprocess_output",
                    excludeTransformParameters=c("FSC-A","FSC-W","FSC-H","Time","Cell_length"))

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






