library(dplyr)
library(tidyr)
library(MetaCyto)
library(ggplot2)


# This example is a light weight example of how to use MetaCyto to  analyze cytometry data from one experiment. 

###################################################
########## Check working directory first ##########
###################################################
wd=getwd()
if(grepl("MetaCyto_Example_One_Experiment$",wd)==F){print("Please change the working directory to 'MetaCyto_Example_One_Experiment'")}


#############################################
########## Collect data Manually  ##########
############################################

# If you want to associate fcs data with some phenotype or genotype data, such as age, gender or SNP, 
# you need to collect the sample information as well, like the one below. 
#  The fcs_files column specifies the location of the fcs files.
# Other two columns specifies the age and gender of each subject corresponding to each fcs file. 
# you can include whatever sample information in the table. 
sample_info=read.csv("Data/sample_info.csv")
print(sample_info)


#############################################
########## Perform preprossessing ##########
############################################
files = list.files("Data",pattern="fcs$",full.names=T)
fcs = preprocessing(fcsFiles=files, 
                    assay = "FCM", 
                    b = 1/150,
                    fileSampleSize=1000, 
                    excludeTransformParameters = c("FSC-A", "FSC-W", "FSC-H", "Time", "Cell_length"))

#######################################
########## Cluster the Data ##########
######################################

# cluster using flowSOM.MC
cluster_list=flowSOM.MC(fcsFrame=fcs,
                        excludeClusterParameters=c("FSC-A","FSC-W","FSC-H","SSC-A","SSC-W","SSC-H","Time"))

# label each clusters
cluster_label=labelCluster(fcsFrame=fcs,
                           clusterList=cluster_list,
                           excludeClusterParameters=c("FSC-A","FSC-W","FSC-H","SSC-A","SSC-W","SSC-H","Time"),
                           minPercent=0.01,
                           labelQuantile=0.95)

# Add a user-defined label (CCR7+ CD8+ T cells) to the cluster_label. This step allows users 
# to include their favorate cell types into the analysis. 
cluster_label$clusterLabel=c(cluster_label$clusterLabel,"CD3+|CD4-|CD8B+|CCR7+")

# Find the hyper-ractangle clusters corresponding to the labels
cluster_list=searchCluster(fcsFrame=fcs,
                      clusterLabel=cluster_label$clusterLabel,
                      cutoff=cluster_label$cutoff)

# derive summary statistics from each clusters
cluster_stats=clusterStats(fcsFrame=fcs, 
                          clusterList=cluster_list$clusterList, 
                          fcsNames=files)


############################################
########## Visualize the results  ##########
############################################

# let's see how the fraction of "CCR7+ CD8+ T cells" change with age
sub_data=subset(cluster_stats,cluster_stats$label=="CD3+|CD4-|CD8B+|CCR7+")
sub_data=inner_join(sub_data,sample_info,by=c("fcs_names"="fcs_files"))
p=ggplot(data=sub_data,aes(x=SUBJECT_AGE,y=fraction))+
  geom_point()+
  geom_smooth(method="lm")+
  ggtitle("CCR7+ CD8+ T cells")
plot(p)

