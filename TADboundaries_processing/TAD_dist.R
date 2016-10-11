# This script intends to display the distribution of TAD, gaps and boundaries length

# Cyril Matthey-Doret 
# 07.10.2016
######################################

setwd("/home/cyril/Documents/First_step/data/")
#setwd("/Users/cmatthe5/Documents/First_step/data/")

# Visualizing the length of RNAs, TADs and boundaries (at 10% threshold)
TADb10 <- read.table("TAD/TAD_boundaries10.bed")
TADb10 <- TADb10[,-4]
TAD <- read.table("TAD/GM12878_TAD_domains.bed")
TAD <- TAD[,-4]
colnames(TADb10) = colnames(TAD)<-c("chr", "start","end","ID")
TADb_lincRNA<-read.table("linc_RNA/TADbound-lincRNA10.bed")
TADb_pcgenes <- read.table("pc_genes/TADbound-pcgene10.bed")
colnames(TADb_lincRNA) = colnames(TADb_pcgenes)<-c("chr", "start","end","ID", "strand")
par(mfrow=c(4,1))
hist(log10(TADb10$end - TADb10$start),xlim=c(2,7))
hist(log10(TAD$end - TAD$start),xlim=c(2,7))
hist(log10(TADb_pcgenes$end - TADb_pcgenes$start),xlim=c(2,7))
hist(log10(TADb_lincRNA$end - TADb_lincRNA$start),xlim=c(2,7))

