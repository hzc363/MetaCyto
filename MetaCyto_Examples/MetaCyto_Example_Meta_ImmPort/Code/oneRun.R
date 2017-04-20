codes=list.files("Code",pattern="step",full.names=T)
sapply(codes,source)