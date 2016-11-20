setwd("/Users/cmatthe5/Documents/First_step/data/")
setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
hicboundfull<-read.table("TAD/hicboundaries/HICBOUND3.txt",header=T)
library(dplyr)

merged<-bind_rows(hicboundfull[,c("chr","Lbound.start","Lbound.end")],hicboundfull[,c("chr","Rbound.start","Rbound.end")])

merged$Lbound.start[is.na(merged$Lbound.start)] <- merged$Rbound.start[is.na(merged$Lbound.start)]
merged$Lbound.end[is.na(merged$Lbound.end)] <- merged$Rbound.end[is.na(merged$Lbound.end)]
merged <- merged[,c(1,2,3)]
options(scipen=999)
write.table(merged,"TAD/hicboundaries/ucsc_chr3.txt",col.names=F,row.names=F,sep="\t",quote=F)
