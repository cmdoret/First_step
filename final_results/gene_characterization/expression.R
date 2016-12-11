# The purpose of this script is to generate a summarized table containing PCG and lincRNAs with their associated median expression levels
# Previously computed from the expression matrix (available on ENCODE) in 4 different cell types, along with their respective overlap 
# status with enhancers and promoters
# Cyril Matthey-Doret, 18.11.2016

library(ggplot2)
library(gridExtra)

######################################################

#Loading data:
setwd("/Users/cmatthe5/Documents/First_step/data/")
#setwd("/home/cyril/Documents/First_step/data/")
setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
bedcol <- c("chr", "start", "end", "gene", "strand")
exp_lincRNA <- read.table("expression/HTseq.count.encode.matrix.RPKM.avg.txt", header = T)
exp_pcgene <- read.table("expression/HTseq.count.encode.matrix.RPKM.avg.txt", header = T)
# Loading lincRNAs sets
ne_linc <- read.table("enhancer_bound/all_combinations/ne_linc_pr.bed")  # Overlap no enhancer marks, does not take promoters into account
e_linc <- read.table("enhancer_bound/all_combinations/e_linc_pr.bed")  # Overlap enhancer marks, does not take promoters into account
ne.p_linc <- read.table("enhancer_bound/all_combinations/ne.p_linc_pr.bed")  # Overlap promoter marks, but no enhancer marks
e.np_linc <- read.table("enhancer_bound/all_combinations/e.np_linc_pr.bed")  # Overlap enhancer marks, but no promoter marks
e.p_linc <- read.table("enhancer_bound/all_combinations/e.p_linc_pr.bed")  # Overlap promoter marks, and enhancer marks
ne.np_linc <- read.table("enhancer_bound/all_combinations/ne.np_linc_pr.bed")  # Overlap neither promoter marks, nor enhancer marks

# Loading pcgenes sets
ne_pc <- read.table("enhancer_bound/all_combinations/ne_pc_pr.bed")  # Overlap no enhancer marks, does not take promoters into account
e_pc <- read.table("enhancer_bound/all_combinations/e_pc_pr.bed")  # Overlap enhancer marks, does not take promoters into account
ne.p_pc <- read.table("enhancer_bound/all_combinations/ne.p_pc_pr.bed")  # Overlap promoter marks, but no enhancer marks
e.np_pc <- read.table("enhancer_bound/all_combinations/e.np_pc_pr.bed")  # Overlap enhancer marks, but no promoter marks
e.p_pc <- read.table("enhancer_bound/all_combinations/e.p_pc_pr.bed")  # Overlap promoter marks, and enhancer marks
ne.np_pc <- read.table("enhancer_bound/all_combinations/ne.np_pc_pr.bed")  # Overlap neither promoter marks, nor enhancer marks


#adding colnames
colnames(ne_linc)= colnames(e_linc)=colnames(ne.p_linc)=colnames(e.np_linc)=colnames(e.p_linc)=
  colnames(ne.np_linc)=colnames(ne_pc)=colnames(e_pc)=colnames(ne.p_pc)=colnames(e.np_pc)=colnames(e.p_pc)=
  colnames(ne.np_pc) <-c("chr", "start", "end", "gene", "strand")
colnames(exp_pcgene)=colnames(exp_lincRNA) <- c("gene", "expression")

#Splitting expression levels data into categories
exp_ne_linc <-exp_lincRNA[exp_lincRNA$gene %in% ne_linc$gene,]
exp_e_linc <-exp_lincRNA[exp_lincRNA$gene %in% e_linc$gene,]
exp_e.p_linc <-exp_lincRNA[exp_lincRNA$gene %in% e.p_linc$gene,]
exp_ne.p_linc <-exp_lincRNA[exp_lincRNA$gene %in% ne.p_linc$gene,]
exp_e.np_linc <-exp_lincRNA[exp_lincRNA$gene %in% e.np_linc$gene,]
exp_ne.np_linc <-exp_lincRNA[exp_lincRNA$gene %in% ne.np_linc$gene,]

exp_ne_pc <-exp_pcgene[exp_pcgene$gene %in% ne_pc$gene,]
exp_e_pc <-exp_pcgene[exp_pcgene$gene %in% e_pc$gene,]
exp_e.p_pc <-exp_pcgene[exp_pcgene$gene %in% e.p_pc$gene,]
exp_ne.p_pc <-exp_pcgene[exp_pcgene$gene %in% ne.p_pc$gene,]
exp_e.np_pc <-exp_pcgene[exp_pcgene$gene %in% e.np_pc$gene,]
exp_ne.np_pc <-exp_pcgene[exp_pcgene$gene %in% ne.np_pc$gene,]


#building large, single dataframe to make data more convenient to use.

whole_exp <- rbind(cbind(exp_e_linc,promoter=rep("both"),enhancer=rep("+"),gentype=rep("linc")), 
                   cbind(exp_ne_linc,promoter=rep("both"),enhancer=rep("-"),gentype=rep("linc")),
                   cbind(exp_e.np_linc,promoter=rep("-"),enhancer=rep("+"),gentype=rep("linc")),
                   cbind(exp_ne.p_linc,promoter=rep("+"),enhancer=rep("-"),gentype=rep("linc")),
                   cbind(exp_e.p_linc,promoter=rep("+"),enhancer=rep("+"),gentype=rep("linc")),
                   cbind(exp_ne.np_linc,promoter=rep("-"),enhancer=rep("-"),gentype=rep("linc")),
                   cbind(exp_e_pc,promoter=rep("both"),enhancer=rep("+"),gentype=rep("pc")), 
                   cbind(exp_ne_pc,promoter=rep("both"),enhancer=rep("-"),gentype=rep("pc")),
                   cbind(exp_e.np_pc,promoter=rep("-"),enhancer=rep("+"),gentype=rep("pc")),
                   cbind(exp_ne.p_pc,promoter=rep("+"),enhancer=rep("-"),gentype=rep("pc")),
                   cbind(exp_e.p_pc,promoter=rep("+"),enhancer=rep("+"),gentype=rep("pc")),
                   cbind(exp_ne.np_pc,promoter=rep("-"),enhancer=rep("-"),gentype=rep("pc")))
write.table(x = whole_exp,file = "expression/enhancer_promoter/whole_exp.txt",quote = F,sep = "\t",row.names = F,col.names = T)

