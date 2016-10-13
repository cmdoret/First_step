
# This script intends to count the number of lincRNAs/protein-coding genes overlapping each bin
# of the merged TADs.
# Cyril
# 13.10.2016
##################################

# Loading data:
setwd("/home/cyril/Documents/First_step/data/")
#setwd("/Users/cmatthe5/Documents/First_step/data/")
TADbins <- read.table("TAD/merged/merged_TADbins.txt")
lincRNA <- read.table("linc_RNA/LCL.expressed.lincRNA.bed")
colnames(lincRNA) <- c("chr","start","end","ID","strand")
colnames(TADbins) <- c("ID","chr","pos","start","end")
over_pc <- read.table("../TAD-lincRNA-overlap/bin-pcgenes_overlap.bed")
over_lincRNA <- read.table("../TAD-lincRNA-overlap/bin-lincRNA_overlap.bed")
#================================

# Counting overlaps:

TADbins <- cbind(TADbins, lincRNA=rep(0,length(TADbins$ID)),pcgene=rep(0,length(TADbins$ID)))

for(id in row(TADbins)){
  print(id)
}