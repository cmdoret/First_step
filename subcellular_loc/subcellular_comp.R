# The purpose of this script is to compare the sucellular localization of TADbound and non-TADbound lincRNAs.
# It also compares the subcellular localizatioin of TADbound and non-TADbound protein coding genes.

# Cyril Matthey-Doret, 05.10.2016

######################################################
library(ggplot2)
library(gridExtra)
#Loading data:
#setwd("/home/cyril/Documents/First_step/data/")
setwd("/home/cyril/Documents/Master/sem_1/First step/data/")
bedcol <- c("chr", "start", "end", "gene", "strand")

loc_lincRNA <- read.table("all.lincRNA.GM12878.subcellular.ratio.txt", header = T)
loc_pcgene <- read.table("all.pcgene.GM12878.subcellular.ratio.txt", header = T)
Tb_lincRNA <- read.table("TADbound-lincRNA.bed")
Tb_pc <- read.table("TADbound-pcgene.bed")
nTb_lincRNA <- read.table("nonTADbound-lincRNA.bed")
nTb_pc <- read.table("nonTADbound-pcgene.bed")
colnames(nTb_pc)=colnames(nTb_lincRNA)=colnames(Tb_lincRNA)=colnames(Tb_pc)<-bedcol #putting informative colnames to bed files

#Splitting localisation data into TADbound (Tb) and non-TADbound (nTb)
loc_Tb_lincRNA <-loc_lincRNA[loc_lincRNA$gene %in% Tb_lincRNA$gene,]
loc_nTb_lincRNA <-loc_lincRNA[loc_lincRNA$gene %in% nTb_lincRNA$gene,]
loc_Tb_pc <-loc_pcgene[loc_pcgene$gene %in% Tb_pc$gene,]
loc_nTb_pc <-loc_pcgene[loc_pcgene$gene %in% nTb_pc$gene,]

#Comparing overall mean ratio between Tb and nTb:
#lincRNAs:
hist_linc <-ggplot()+
  geom_histogram(data=loc_Tb_lincRNA, aes(x=ratio, y=..density..), fill="#0000dd", alpha=0.5, bins = 300)+
  geom_histogram(data=loc_nTb_lincRNA, aes(x=ratio, y=..density..), fill="#bb0000", alpha=0.5, bins = 300)+
  ggtitle("lincRNAs")
summary(loc_Tb_lincRNA$ratio);summary(loc_nTb_lincRNA$ratio)

#Protein-coding genes:
hist_pc<-ggplot()+
  geom_histogram(data=loc_Tb_pc, aes(x=ratio, y=..density..), fill="#0000dd", alpha=0.5, bins = 300)+
  geom_histogram(data=loc_nTb_pc, aes(x=ratio, y=..density..), fill="#bb0000", alpha=0.5, bins = 300)+
  ggtitle("Protein-coding genes")
summary(loc_Tb_lincRNA$ratio);summary(loc_nTb_lincRNA$ratio)

grid.arrange(hist_linc, hist_pc)
