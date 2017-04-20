labelUnifier=function(clusterLabel){
  clusterLabel=strsplit(as.character(clusterLabel),split="\\|")
  clusterLabel=lapply(clusterLabel,sort)
  clusterLabel=sapply(clusterLabel,paste,collapse="|")
  clusterLabel=toupper(clusterLabel)
  return(clusterLabel)
}
labelCombiner=function(labels){
  label1=labels[1]
  label2=labels[2]
  if(nchar(label1)!=nchar(label2)){return(NA)}
  m1=gsub("\\+|-","",label1)
  m2=gsub("\\+|-","",label2)
  if(m1!=m2){return(NA)}
  label1=strsplit(label1,"\\|")[[1]]
  label2=strsplit(label2,"\\|")[[1]]
  s1=sapply(label1, function(x){substring(x,nchar(x))})
  s2=sapply(label2, function(x){substring(x,nchar(x))})

  if(sum(s1!=s2)!=1){return(NA)}
  w=which(s1==s2)
  label=label1[w]
  label=paste(label,collapse="|")
  return(label)
}
HC.assign=function(HC,minimumClusterSizePercent){
  CL=list()
  j=0
  NL=min(100,nrow(HC$merge))
  for(k in 2:100){
    t1=cutree(HC,k)
    w=lapply(unique(t1),function(x){which(t1==x)})
    w=w[sapply(w, length)>minimumClusterSizePercent*(nrow(HC$merge)+1)]

    if(length(w)==0){break}
    CL_new=c(CL,w)
    if(length(unique(CL_new))==length(unique(CL))){j=j+1}else{j=0}
    CL=unique(CL_new)
    if(j>20){break}
  }
  return(CL)
}
baselineCut=function(x){
  breaks=seq(min(x),max(x),length.out=100)
  h=hist(x,breaks,plot=FALSE)
  #h2=h$counts[h$counts>1]
  baseline=max(4,mean(h$counts)/10)
  bound=range(which(h$counts>baseline))
  x=x[x>breaks[bound[1]]&x<breaks[bound[2]]]
  return(x)
}
trisect=function(x){
  if(length(x)>2000){x=sample(x,2000)}
  x=baselineCut(x)
  D=dist(x)
  cutoff1=findCutoff(x)

  t1=findCutoff(x[x>cutoff1],useBL=FALSE,minX=NULL)
  cluster=.bincode(x,sort(c(-Inf,cutoff1,t1,+Inf),decreasing=FALSE))
  ss <- cluster::silhouette(cluster,D)
  s1=mean(ss[, 3])

  t2=findCutoff(x[x<cutoff1],useBL=FALSE,minX=NULL)
  cluster=.bincode(x,sort(c(-Inf,cutoff1,t2,+Inf),decreasing=FALSE))
  ss <- cluster::silhouette(cluster,D)
  s2=mean(ss[, 3])

  if(s2>s1){cutoff2=t2}else{cutoff2=t1}
  return(sort(c(cutoff1,cutoff2),decreasing=TRUE))
}
