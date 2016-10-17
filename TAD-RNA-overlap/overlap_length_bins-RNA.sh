#!bin/bash

# This script intends to measure all overlaps between TAD bins and lincRNAs or protein coding genes. 
# It runs with 3 different binwidth and there is no requirement for the minimum overlap length.
# Cyril Matthey-Doret
# 13.10.2016


for bw in 1 5 10;
do
bedtools intersect -a ../data/TAD/merged/short_TADbins"$bw".bed -b ../data/pc_genes/LCL.expressed.pcgene.bed -wao > bin-pcgenes_overlap"$bw"_length.bed;
bedtools intersect -a ../data/TAD/merged/short_TADbins"$bw".bed -b ../data/linc_RNA/LCL.expressed.lincRNA.bed -wao > bin-lincRNA_overlap"$bw"_length.bed;
done

