# The purpose of this script is to compare the sucellular localization of TADbound and non-TADbound lincRNAs.
# It also compares the subcellular localizatioin of TADbound and non-TADbound protein coding genes.

# Cyril Matthey-Doret, 05.10.2016

######################################################

#Loading data:
setwd("/home/cyril/Documents/First_step/data/")
loc_lincRNA <- read.table("all.lincRNA.GM12878.subcellular.ratio.txt", header = T)
loc_pcgene <- read.table("all.pcgene.GM12878.subcellular.ratio.txt", header = T)
Tb_lincRNA <- read.table("TADbound-lincRNA.bed")
Tb_pc <- read.table("TADbound-pcgene.bed")
nTb_lincRNA <- read.table("nonTADbound-lincRNA.bed")
nTb_pc <- read.table("nonTADbound-pcgene.bed")

#Splitting localisation data into TADbound (Tb) and non-TADbound (nTb)

