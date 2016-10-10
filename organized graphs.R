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
exp_lincRNA <- read.table("expression/LCL.lincRNA.expression.txt", header = F)
exp_pcgene <- read.table("expression/LCL.pcgene.expression.txt", header = F)
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

#adding colnames
colnames(Tb_lincRNA5)= colnames(Tb_lincRNA10)=colnames(Tb_lincRNA20)=colnames(nTb_lincRNA5)=colnames(nTb_lincRNA10)=
  colnames(nTb_lincRNA20)=colnames(Tb_pc5)=colnames(Tb_pc10)=colnames(Tb_pc20)=colnames(nTb_pc5)=colnames(nTb_pc10)=
  colnames(nTb_pc20) <-c("chr", "start", "end", "gene", "strand")
colnames(exp_pcgene)=colnames(exp_lincRNA) <- c("gene", "expression")

#Splitting expression levels data into TADbound (Tb) and non-TADbound (nTb)
exp_Tb_lincRNA5 <-exp_lincRNA[exp_lincRNA$gene %in% Tb_lincRNA5$gene,]
exp_Tb_lincRNA10 <-exp_lincRNA[exp_lincRNA$gene %in% Tb_lincRNA10$gene,]
exp_Tb_lincRNA20 <-exp_lincRNA[exp_lincRNA$gene %in% Tb_lincRNA20$gene,]

exp_nTb_lincRNA5 <-exp_lincRNA[exp_lincRNA$gene %in% nTb_lincRNA5$gene,]
exp_nTb_lincRNA10 <-exp_lincRNA[exp_lincRNA$gene %in% nTb_lincRNA10$gene,]
exp_nTb_lincRNA20 <-exp_lincRNA[exp_lincRNA$gene %in% nTb_lincRNA20$gene,]

exp_Tb_pc5 <-exp_pcgene[exp_pcgene$gene %in% Tb_pc5$gene,]
exp_Tb_pc10 <-exp_pcgene[exp_pcgene$gene %in% Tb_pc10$gene,]
exp_Tb_pc20 <-exp_pcgene[exp_pcgene$gene %in% Tb_pc20$gene,]

exp_nTb_pc5 <-exp_pcgene[exp_pcgene$gene %in% nTb_pc5$gene,]
exp_nTb_pc10 <-exp_pcgene[exp_pcgene$gene %in% nTb_pc10$gene,]
exp_nTb_pc20 <-exp_pcgene[exp_pcgene$gene %in% nTb_pc20$gene,]

#Writing into bed files for further use.
#write.table(exp_Tb_lincRNA5,file = "expression/exp_Tb_lincRNA5.txt",sep="\t",quote = F,col.names = F,row.names = F)
#write.table(exp_Tb_lincRNA10,file = "expression/exp_Tb_lincRNA10.txt",sep="\t",quote = F,col.names = F,row.names = F)
#write.table(exp_Tb_lincRNA20,file = "expression/exp_Tb_lincRNA20.txt",sep="\t",quote = F,col.names = F,row.names = F)

#write.table(exp_nTb_lincRNA5,file = "expression/exp_nTb_lincRNA5.txt",sep="\t",quote = F,col.names = F,row.names = F)
#write.table(exp_nTb_lincRNA10,file = "expression/exp_nTb_lincRNA10.txt",sep="\t",quote = F,col.names = F,row.names = F)
#write.table(exp_nTb_lincRNA20,file = "expression/exp_nTb_lincRNA20.txt",sep="\t",quote = F,col.names = F,row.names = F)

#write.table(exp_Tb_pc5,file = "expression/exp_Tb_pc5.txt",sep="\t",quote = F,col.names = F,row.names = F)
#write.table(exp_Tb_pc10,file = "expression/exp_Tb_pc10.txt",sep="\t",quote = F,col.names = F,row.names = F)
#write.table(exp_Tb_pc20,file = "expression/exp_Tb_pc20.txt",sep="\t",quote = F,col.names = F,row.names = F)

#write.table(exp_nTb_pc5,file = "expression/exp_nTb_pc5.txt",sep="\t",quote = F,col.names = F,row.names = F)
#write.table(exp_nTb_pc10,file = "expression/exp_nTb_pc10.txt",sep="\t",quote = F,col.names = F,row.names = F)
#write.table(exp_nTb_pc20,file = "expression/exp_nTb_pc20.txt",sep="\t",quote = F,col.names = F,row.names = F)

#building large, single dataframe to make data more convenient.

CData <-function(df){
  return(cbind(df,fact=rep(deparse(substitute(df)),length(df[,1]))))
}

whole_exp <- rbind(cbind(exp_Tb_pc5,threshold=rep("5"),TAD=rep("Tb"),gentype=rep("pc")), 
                   cbind(exp_Tb_pc10,threshold=rep("10"),TAD=rep("Tb"),gentype=rep("pc")), 
                   cbind(exp_Tb_pc20,threshold=rep("20"),TAD=rep("Tb"),gentype=rep("pc")), 
                   cbind(exp_nTb_pc5,threshold=rep("5"),TAD=rep("nTb"),gentype=rep("pc")), 
                   cbind(exp_nTb_pc10,threshold=rep("10"),TAD=rep("nTb"),gentype=rep("pc")), 
                   cbind(exp_nTb_pc20,threshold=rep("20"),TAD=rep("nTb"),gentype=rep("pc")),
                   cbind(exp_Tb_lincRNA5,threshold=rep("5"),TAD=rep("Tb"),gentype=rep("lincRNA")), 
                   cbind(exp_Tb_lincRNA10,threshold=rep("10"),TAD=rep("Tb"),gentype=rep("lincRNA")), 
                   cbind(exp_Tb_lincRNA20,threshold=rep("20"),TAD=rep("Tb"),gentype=rep("lincRNA")), 
                   cbind(exp_nTb_lincRNA5,threshold=rep("5"),TAD=rep("nTb"),gentype=rep("lincRNA")), 
                   cbind(exp_nTb_lincRNA10,threshold=rep("10"),TAD=rep("nTb"),gentype=rep("lincRNA")), 
                   cbind(exp_nTb_lincRNA20,threshold=rep("20"),TAD=rep("nTb"),gentype=rep("lincRNA")))
#====================================================

#Visualizing data:
library(plyr)
short_med<-function(x){return(round(median(log10(x)),3))}
short_wilcox <- function(x,y){return(format(wilcox.test(x, y)$p.value,digits=3))}

#dataf frame for annotation of median
med.fac = ddply(whole_exp, .(gentype, threshold,TAD), function(.d)
  data.frame(x=median(round(log10(.d$expression),3))))
#data frame for annotation of p-value
wil.fac = ddply(whole_exp, .(gentype, threshold,TAD), function(.d)
  data.frame(x=short_wilcox(.d$expression~.d$TAD)))


ggplot(data=whole_exp)+
  facet_grid(gentype~threshold)+
  geom_boxplot(aes(x=TAD, y=log10(expression),fill=TAD),notch=T)+
  geom_text(data=med.fac, aes(x=TAD, y=-1, label=x), 
            colour="black", inherit.aes=FALSE, parse=FALSE)+
  geom_text(data=wil.fac, aes(x=1.5, y=1, label=x), 
            colour="black", inherit.aes=FALSE, parse=FALSE)+
  theme_bw()

#====================================================

#Summary statistics:

summary(); summary()
summary(exp_Tb_pc$expression); summary(exp_nTb_pc$expression)
