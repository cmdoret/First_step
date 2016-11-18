# The purpose of this script is to compare the sucellular localization of TADbound and non-TADbound lincRNAs.
# It also compares the subcellular localizatioin of TADbound and non-TADbound protein coding genes.

# Cyril Matthey-Doret, 05.10.2016

library(ggplot2)
library(gridExtra)

######################################################

#Loading data:
setwd("/Users/cmatthe5/Documents/First_step/data/")
#setwd("/home/cyril/Documents/First_step/data/")
#setwd("/home/cyril/Documents/Master/sem_1/First step/data/")
bedcol <- c("chr", "start", "end", "gene", "strand")
loc_lincRNA <- read.table("subcell_loc/all.lincRNA.GM12878.subcellular.ratio.txt", header = T)
loc_pcgene <- read.table("subcell_loc/all.pcgene.GM12878.subcellular.ratio.txt", header = T)
#loading genes sets
elinc <- read.table("enhancer_bound/elinc_prb.bed")
plinc<- read.table("enhancer_bound/plinc_prb.bed")
epc<- read.table("enhancer_bound/epc_prb.bed")
ppc<- read.table("enhancer_bound/ppc_prb.bed")

colnames(elinc)=colnames(plinc)=colnames(epc)=
  colnames(ppc) <-c("chr", "start", "end", "gene", "strand")


#Splitting localisation data into TADbound (Tb) and non-TADbound (nTb)

loc_elinc <-loc_lincRNA[loc_lincRNA$gene %in% elinc$gene,]
loc_plinc <-loc_lincRNA[loc_lincRNA$gene %in% plinc$gene,]
loc_epc <-loc_pcgene[loc_pcgene$gene %in% epc$gene,]
loc_ppc <-loc_pcgene[loc_pcgene$gene %in% ppc$gene,]

#============================================================

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
#It doesn't seem TADbound lincRNAs/proteins are more localized in the nucleus/cytoplasm than non-TADbound ones.
boxplot(notch=T, log10(loc_Tb_lincRNA$ratio), log10(loc_nTb_lincRNA$ratio))
wilcox.test(loc_Tb_lincRNA$ratio, loc_nTb_lincRNA$ratio)
#==========================================================

#building large, single dataframe to make data more convenient.

whole_loc <- rbind(cbind(loc_elinc,assoc=rep("e"),gentype=rep("lincRNA")), 
                   cbind(loc_plinc,assoc=rep("p"),gentype=rep("lincRNA")), 
                   cbind(loc_epc,assoc=rep("e"),gentype=rep("pc")), 
                   cbind(loc_ppc,assoc=rep("p"),gentype=rep("pc")))
write.table(x = whole_loc,file = "subcell_loc/enhancer_bound/whole_loc.txt",quote = F,sep = "\t",row.names = F,col.names = T)


