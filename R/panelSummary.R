#' Summarize markers in panels.
#'
#' A function that summarizes markers in cytometry panels.
#' @param panelInfo A data frame returned by the collectData function. It should contain all the information outputted by the preprocessing.batch function.
#' @param folder The directory where the output should be written.
#' @param cluster True or False. Used to indicate if the markers and panels should be clustered in the plot.
#' @param plotImage True or False. Used to indicate if a plot summarizing markers in panels should be produced.
#' @param width Used to specify the width of the plot
#' @param height Used to specify the height of the plot
#' @return A dataframe describing what markers are in each panel.
#' @examples
#' fn=system.file("extdata","",package="MetaCyto")
#' fn=list.files(fn,pattern="processed_sample",full.names=TRUE)
#' panel_info=collectData(fn,longform=FALSE)
#' dir.create("Example_Result")
#' PS=panelSummary(panel_info,"Example_Result",cluster=FALSE,width=30,height=20)
#' @export
panelSummary=function(panelInfo,folder,cluster=TRUE,plotImage=TRUE,width=20,height=20){
  panelInfo=unique(panelInfo[,c("study_id","antibodies")])
  ab_list=NULL
  for(i in 1:nrow(panelInfo)){
    t1=strsplit(panelInfo$antibodies[i],split="\\|")[[1]]
    t1=t1[t1!="NA"&is.na(t1)==FALSE]
    if(length(t1)>0){
      ab_list=rbind(ab_list, data.frame("study_id"=panelInfo$study_id[i],"antibodies"=t1,"value"=1))
    }
  }
  ab_list=unique(ab_list)
  #ab_list=unique(ab_list)
  ab_table=tidyr::spread(ab_list,key=antibodies,value=value,fill=0)
  rownames(ab_table)=ab_table$study_id;ab_table=ab_table[,-1]
  ab_table=t(data.matrix(ab_table))
  if(cluster==TRUE){
    d=as.dist(1-cor(ab_table))
    c1=hclust(d)$order
    d=as.dist(1-cor(t(ab_table)))
    r1=hclust(d)$order
    ab_table=ab_table[r1,c1]
  }
  write.csv(ab_table,paste0(folder,"/panel_summary.csv"))

  if(plotImage==TRUE){
    pdf(paste0(folder,"/panel_summary.pdf"),width=width,height=height)
    op <- par(mar = c(30,30,30,10))
    image(z = ab_table, col = c("white","red"), axes = FALSE)
    t1=1/(nrow(ab_table)-1)
    axis(side = 3, labels = rownames(ab_table),las=2,cex.axis=3,
         at = seq(0, 1,by = t1))
    t1=1/(ncol(ab_table)-1)
    axis(side = 2, labels = colnames(ab_table),las=2,cex.axis=3,
         at = seq(0, 1,by = t1))
    box()
    par(op)
    dev.off()
  }
  return(ab_table)
}
