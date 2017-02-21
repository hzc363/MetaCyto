library(MetaCyto)

#################################################
########## Prepare for preprossessing ##########
################################################

# Read the fcs_info file from step 1
fcs_info=read.csv("Data/fcs_info.csv",stringsAsFactors=F,check.names=F)

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


