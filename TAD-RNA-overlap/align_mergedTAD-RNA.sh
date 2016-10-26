#!/bin/bash

#This script intends to align lincRNA and protein coding genes with TAD boundaries of merged TADs. It outputs RNAs which have at least
#25% of their sequence overlapping a TAD boundary.
# Cyril Matthey-Doret, 04.10.2016

for t in 5;
do
bedtools intersect -a ../data/linc_RNA/LCL.expressed.lincRNA.bed -b ../data/TAD/merged/flexibound"$t".bed -f 0.51 -wa -wb > NEWlincRNA_"$t"overlap_flexible_TADb.bed;

bedtools intersect -a ../data/pc_genes/LCL.expressed.pcgene.bed -b ../data/TAD/merged/flexibound"$t".bed -f 0.11 -wa -wb > NEWpcgene_"$t"overlap_flexible_TADb.bed;

done
