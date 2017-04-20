library(MetaCyto)

###############################################
########## excludeClusterParameters  ##########
###############################################

#define parameters that we don't want to cluster
excludeClusterParameters=c("FSC-A","FSC-W","FSC-H","SSC-A","SSC-W","SSC-H","Time",
                           "CELL_LENGTH","DEAD","DNA1","DNA2")


###############################################
########## Cluster the cytometry data  ########
###############################################

# Find and label clusters in the data. The default cluster functions 
# are FlowSOM.MC. Here, flowHC, a hiarachical clustering function is used instead
# to show how non-default functions can be used. 
cluster_label=autoCluster.batch(preprocessOutputFolder="Result/preprocess_output",
                                 excludeClusterParameters=excludeClusterParameters,
                                 labelQuantile=0.95,
                                 clusterFunction=flowHC)

# Add a user-defined label (CCR7+ CD8+ T cells) to the cluster_label. This step allows users 
# to include their favorate cell type into the analysis. 
cluster_label=c(cluster_label,"CD3+|CD4-|CD8+|CCR7+")

# Derive summary statistics for the clusters
# Result will be written out the the directory speficied by the "outpath" argument
searchCluster.batch(preprocessOutputFolder="Result/preprocess_output",
              outpath="Result/search_output",
              clusterLabel=cluster_label)

