# The purpose of this script is to compare the sucellular localization of TADbound and non-TADbound lincRNAs.
# It also compares the subcellular localizatioin of TADbound and non-TADbound protein coding genes.

# Cyril Matthey-Doret, 05.10.2016

library(ggplot2)
library(gridExtra)

######################################################

#Loading data:
setwd("/home/cyril/Documents/First_step/data/")
setwd("/Users/cmatthe5/Documents/First_step/data/")
#setwd("/home/cyril/Documents/Master/sem_1/First step/data/")
bedcol <- c("chr", "start", "end", "gene", "strand")
loc_lincRNA <- read.table("subcell_loc/all.lincRNA.GM12878.subcellular.ratio.txt", header = T)
loc_pcgene <- read.table("subcell_loc/all.pcgene.GM12878.subcellular.ratio.txt", header = T)
#loading lincRNAs sets
ne_linc <- read.table("enhancer_bound/all_combinations/ne_linc_prb.bed")  # Overlap no enhancer marks, does not take promoters into account
e_linc <- read.table("enhancer_bound/all_combinations/e_linc_prb.bed")  # Overlap enhancer marks, does not take promoters into account
ne.p_linc <- read.table("enhancer_bound/all_combinations/ne.p_linc_prb.bed")  # Overlap promoter marks, but no enhancer marks
e.np_linc <- read.table("enhancer_bound/all_combinations/e.np_linc_prb.bed")  # Overlap enhancer marks, but no promoter marks
e.p_linc <- read.table("enhancer_bound/all_combinations/e.p_linc_prb.bed")  # Overlap promoter marks, and enhancer marks
ne.np_linc <- read.table("enhancer_bound/all_combinations/ne.np_linc_prb.bed")  # Overlap neither promoter marks, nor enhancer marks

# Loading pcgenes sets
ne_pc <- read.table("enhancer_bound/all_combinations/ne_pc_prb.bed")  # Overlap no enhancer marks, does not take promoters into account
e_pc <- read.table("enhancer_bound/all_combinations/e_pc_prb.bed")  # Overlap enhancer marks, does not take promoters into account
ne.p_pc <- read.table("enhancer_bound/all_combinations/ne.p_pc_prb.bed")  # Overlap promoter marks, but no enhancer marks
e.np_pc <- read.table("enhancer_bound/all_combinations/e.np_pc_prb.bed")  # Overlap enhancer marks, but no promoter marks
e.p_pc <- read.table("enhancer_bound/all_combinations/e.p_pc_prb.bed")  # Overlap promoter marks, and enhancer marks
ne.np_pc <- read.table("enhancer_bound/all_combinations/ne.np_pc_prb.bed")  # Overlap neither promoter marks, nor enhancer marks


colnames(ne_linc)= colnames(e_linc)=colnames(ne.p_linc)=colnames(e.np_linc)=colnames(e.p_linc)=
  colnames(ne.np_linc)=colnames(ne_pc)=colnames(e_pc)=colnames(ne.p_pc)=colnames(e.np_pc)=colnames(e.p_pc)=
  colnames(ne.np_pc) <-c("chr", "start", "end", "gene", "strand")


#Splitting localisation data into categories

loc_ne_linc <-loc_lincRNA[loc_lincRNA$gene %in% ne_linc$gene,]
loc_e_linc <-loc_lincRNA[loc_lincRNA$gene %in% e_linc$gene,]
loc_e.p_linc <-loc_lincRNA[loc_lincRNA$gene %in% e.p_linc$gene,]
loc_ne.p_linc <-loc_lincRNA[loc_lincRNA$gene %in% ne.p_linc$gene,]
loc_e.np_linc <-loc_lincRNA[loc_lincRNA$gene %in% e.np_linc$gene,]
loc_ne.np_linc <-loc_lincRNA[loc_lincRNA$gene %in% ne.np_linc$gene,]

loc_ne_pc <-loc_pcgene[loc_pcgene$gene %in% ne_pc$gene,]
loc_e_pc <-loc_pcgene[loc_pcgene$gene %in% e_pc$gene,]
loc_e.p_pc <-loc_pcgene[loc_pcgene$gene %in% e.p_pc$gene,]
loc_ne.p_pc <-loc_pcgene[loc_pcgene$gene %in% ne.p_pc$gene,]
loc_e.np_pc <-loc_pcgene[loc_pcgene$gene %in% e.np_pc$gene,]
loc_ne.np_pc <-loc_pcgene[loc_pcgene$gene %in% ne.np_pc$gene,]

#==========================================================

#Building large, single dataframe to make data more convenient.

whole_loc <- rbind(cbind(loc_e_linc,promoter=rep("both"),enhancer=rep("+"),gentype=rep("linc")), 
                   cbind(loc_ne_linc,promoter=rep("both"),enhancer=rep("-"),gentype=rep("linc")),
                   cbind(loc_e.np_linc,promoter=rep("-"),enhancer=rep("+"),gentype=rep("linc")),
                   cbind(loc_ne.p_linc,promoter=rep("+"),enhancer=rep("-"),gentype=rep("linc")),
                   cbind(loc_e.p_linc,promoter=rep("+"),enhancer=rep("+"),gentype=rep("linc")),
                   cbind(loc_ne.np_linc,promoter=rep("-"),enhancer=rep("-"),gentype=rep("linc")),
                   cbind(loc_e_pc,promoter=rep("both"),enhancer=rep("+"),gentype=rep("pc")), 
                   cbind(loc_ne_pc,promoter=rep("both"),enhancer=rep("-"),gentype=rep("pc")),
                   cbind(loc_e.np_pc,promoter=rep("-"),enhancer=rep("+"),gentype=rep("pc")),
                   cbind(loc_ne.p_pc,promoter=rep("+"),enhancer=rep("-"),gentype=rep("pc")),
                   cbind(loc_e.p_pc,promoter=rep("+"),enhancer=rep("+"),gentype=rep("pc")),
                   cbind(loc_ne.np_pc,promoter=rep("-"),enhancer=rep("-"),gentype=rep("pc")))
write.table(x = whole_loc,file = "subcell_loc/enhancer_promoter/whole_loc.txt",quote = F,sep = "\t",row.names = F,col.names = T)
