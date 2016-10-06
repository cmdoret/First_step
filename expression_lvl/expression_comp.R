# The purpose of this script is to compare the median expression level between TADbound and non-TADbound lincRNAs.
# It also compares expression levels between TADbound and non-TADbound protein-coding genes.

# Cyril Matthey-Doret, 05.10.2016

library(ggplot2)
library(gridExtra)

######################################################

#Loading data:
setwd("/home/cyril/Documents/First_step/data/")
#setwd("/home/cyril/Documents/Master/sem_1/First step/data/")
bedcol <- c("chr", "start", "end", "gene", "strand")

exp_lincRNA <- read.table("LCL.lincRNA.expression.txt", header = F)
exp_pcgene <- read.table("LCL.pcgene.expression.txt", header = F)
Tb_lincRNA <- read.table("TADbound-lincRNA.bed")
Tb_pc <- read.table("TADbound-pcgene.bed")
nTb_lincRNA <- read.table("nonTADbound-lincRNA.bed")
nTb_pc <- read.table("nonTADbound-pcgene.bed")
colnames(nTb_pc)=colnames(nTb_lincRNA)=colnames(Tb_lincRNA)=colnames(Tb_pc)<-bedcol #putting informative colnames to bed files
colnames(exp_pcgene)=colnames(exp_lincRNA) <- c("gene", "expression")

#Splitting expression levels data into TADbound (Tb) and non-TADbound (nTb)
exp_Tb_lincRNA <-exp_lincRNA[exp_lincRNA$gene %in% Tb_lincRNA$gene,]
exp_nTb_lincRNA <-exp_lincRNA[exp_lincRNA$gene %in% nTb_lincRNA$gene,]
exp_Tb_pc <-exp_pcgene[exp_pcgene$gene %in% Tb_pc$gene,]
exp_nTb_pc <-exp_pcgene[exp_pcgene$gene %in% nTb_pc$gene,]

#====================================================

#Visualizing data:

hist_linc <-ggplot()+
geom_histogram(data=exp_Tb_lincRNA, aes(x=log10(expression), y=..density..), fill="#0000dd", alpha=0.5, bins = 60)+
  geom_histogram(data=exp_nTb_lincRNA, aes(x=log10(expression), y=..density..), fill="#bb0000", alpha=0.5, bins = 60)+
  ggtitle("lincRNAs")
hist_pc <-ggplot()+
  geom_histogram(data=exp_Tb_pc, aes(x=log10(expression), y=..density..), fill="#0000dd", alpha=0.5, bins = 60)+
  geom_histogram(data=exp_nTb_pc, aes(x=log10(expression), y=..density..), fill="#bb0000", alpha=0.5, bins = 60)+
  ggtitle("lincRNAs")
grid.arrange(hist_linc, hist_pc)
#boxplots
boxplot(log10(exp_nTb_lincRNA$expression), log10(exp_Tb_pc$expression),notch = T)
wilcox.test(exp_Tb_lincRNA$expression, exp_nTb_lincRNA$expression  )
#====================================================

#Summary statistics:

summary(); summary()
summary(exp_Tb_pc$expression); summary(exp_nTb_pc$expression)
