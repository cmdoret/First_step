# The purpose of this script is to compare the median expression level between TADbound and non-TADbound lincRNAs.
# It also compares expression levels between TADbound and non-TADbound protein-coding genes.

# Cyril Matthey-Doret, 05.10.2016

library(ggplot2)
library(gridExtra)

######################################################

#Loading data:
#setwd("/home/cyril/Documents/First_step/data/")
#setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
setwd("/Users/cmatthe5/Documents/First_step/data/")
bedcol <- c("chr", "start", "end", "gene", "strand")
exp_lincRNA <- read.table("expression/LCL.lincRNA.expression.txt", header = F)
exp_pcgene <- read.table("expression/LCL.pcgene.expression.txt", header = F)
test <- read.table("pc_genes/LCL.expressed.pcgene.bed")

#loading enhancer-bound and promoter bound lincRNAs and protein-coding genes
elinc_prb <- read.table("enhancer_bound/elinc_prb.bed")
plinc_prb <- read.table("enhancer_bound/plinc_prb.bed")
epc_prb <- read.table("enhancer_bound/epc_prb.bed")
ppc_prb <- read.table("enhancer_bound/ppc_prb.bed")

#adding colnames
colnames(plinc_prb)=colnames(elinc_prb)=colnames(epc_prb)=
  colnames(ppc_prb) <- bedcol
colnames(exp_pcgene)=colnames(exp_lincRNA) <- c("gene", "expression")

#Splitting expression levels

exp_elinc <- exp_lincRNA[exp_lincRNA$gene %in% elinc_prb$gene,]
exp_plinc <- exp_lincRNA[exp_lincRNA$gene %in% plinc_prb$gene,]
exp_epc <- exp_pcgene[exp_pcgene$gene %in% epc_prb$gene,]
exp_ppc <-exp_pcgene[exp_pcgene$gene %in% ppc_prb$gene,]


#building large, single dataframe to make data more convenient.

whole_exp <- rbind(cbind(exp_elinc,assoc=rep("e"),gentype=rep("lincRNA")), 
                   cbind(exp_plinc,assoc=rep("p"),gentype=rep("lincRNA")), 
                   cbind(exp_epc,assoc=rep("e"),gentype=rep("pc")), 
                   cbind(exp_ppc,assoc=rep("p"),gentype=rep("pc")))
write.table(x = whole_exp,file = "expression/enhancer_bound/whole_exp.txt",quote = F,sep = "\t",row.names = F,col.names = T)


