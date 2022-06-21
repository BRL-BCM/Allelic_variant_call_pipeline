# setwd("../bedgraphs/unionbedg")
args = commandArgs(trailingOnly=TRUE)
STUDY=as.character(args[1])
# STUDY="EXR-TTUSC1gCrGDH-AN"

bedg.bySample=read.table(file = sprintf("%s.bed",STUDY),header = T)
bedg.byStudy=rowSums(bedg.bySample[,4:ncol(bedg.bySample)])
bedg.byStudy=cbind(bedg.bySample[,1:3],bedg.byStudy)
colnames(bedg.byStudy)[colnames(bedg.byStudy)=="bedg.byStudy"]=STUDY
write.table(bedg.byStudy,
            file = sprintf("%s.sum.bed",STUDY),
            quote = F,sep = "\t",
            row.names = F, col.names = T)
