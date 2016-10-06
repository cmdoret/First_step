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
#loading TAD-bound lincRNAs sets
Tb_lincRNA5 <- read.table("linc_RNA/TADbound-lincRNA5.bed")
Tb_lincRNA10 <- read.table("linc_RNA/TADbound-lincRNA10.bed")
Tb_lincRNA20 <- read.table("linc_RNA/TADbound-lincRNA20.bed")
#loading TAD-bound pcgenes sets
Tb_pc5 <- read.table("pc_genes/TADbound-pcgene5.bed")
Tb_pc10 <- read.table("pc_genes/TADbound-pcgene10.bed")
Tb_pc20 <- read.table("pc_genes/TADbound-pcgene20.bed")
#loading non-TAD-bound lincRNAs sets
nTb_lincRNA5 <- read.table("linc_RNA/nonTADbound-lincRNA5.bed")
nTb_lincRNA10 <- read.table("linc_RNA/nonTADbound-lincRNA10.bed")
nTb_lincRNA20 <- read.table("linc_RNA/nonTADbound-lincRNA20.bed")
#loading non-TAD-bound pcgenes sets
nTb_pc5 <- read.table("pc_genes/nonTADbound-pcgene5.bed")
nTb_pc10 <- read.table("pc_genes/nonTADbound-pcgene10.bed")
nTb_pc20 <- read.table("pc_genes/nonTADbound-pcgene20.bed")

ChangeNames <- function(x) {
  names(x) <- c("chr", "start", "end", "gene", "strand")
  return(x)
}
dfs<-list(Tb_lincRNA5, Tb_lincRNA10, Tb_lincRNA20, nTb_lincRNA5, nTb_lincRNA10, nTb_lincRNA20, Tb_pc5, Tb_pc10, Tb_pc20, nTb_pc5, nTb_pc10, nTb_pc20)
dfs <- lapply(dfs, ChangeNames)

for(i in 1:length(x)){
  colnames(x[i])
}
colnames(nTb_pc)=colnames(nTb_lincRNA)=colnames(Tb_lincRNA)=colnames(Tb_pc)<-bedcol #putting informative colnames to bed files
colnames(exp_pcgene)=colnames(exp_lincRNA) <- c("gene", "expression")

#Splitting expression levels data into TADbound (Tb) and non-TADbound (nTb)
exp_Tb_lincRNA <-exp_lincRNA[exp_lincRNA$gene %in% Tb_lincRNA$gene,]
exp_nTb_lincRNA <-exp_lincRNA[exp_lincRNA$gene %in% nTb_lincRNA$gene,]
exp_Tb_pc <-exp_pcgene[exp_pcgene$gene %in% Tb_pc$gene,]
exp_nTb_pc <-exp_pcgene[exp_pcgene$gene %in% nTb_pc$gene,]

#Writing into bed files for further use.
write.table(exp_Tb_lincRNA5,file = "exp_Tb_lincRNA5.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(exp_Tb_lincRNA10,file = "exp_Tb_lincRNA10.txt",sep="\t",quote = F,col.names = F,row.names = F)
vwrite.table(exp_Tb_lincRNA20,file = "exp_Tb_lincRNA20.txt",sep="\t",quote = F,col.names = F,row.names = F)

write.table(exp_nTb_lincRNA5,file = "exp_nTb_lincRNA5.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(exp_nTb_lincRNA10,file = "exp_nTb_lincRNA10.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(exp_nTb_lincRNA20,file = "exp_nTb_lincRNA20.txt",sep="\t",quote = F,col.names = F,row.names = F)

write.table(exp_Tb_pc5,file = "exp_Tb_pc5.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(exp_Tb_pc10,file = "exp_Tb_pc10.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(exp_Tb_pc20,file = "exp_Tb_pc20.txt",sep="\t",quote = F,col.names = F,row.names = F)

write.table(exp_nTb_pc5,file = "exp_nTb_pc5.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(exp_nTb_pc10,file = "exp_nTb_pc10.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(exp_nTb_pc20,file = "exp_nTb_pc20.txt",sep="\t",quote = F,col.names = F,row.names = F)

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
