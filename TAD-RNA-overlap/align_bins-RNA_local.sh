#!bin/bash

# This scripts intend to count the number of overlaps between each TAD bin 
# and at least 25% of the transcripts of protein-coding genes/lincRNAs.
# This is the improved version of the script that requires bedtools v.2.25 (vital-IT only has 2.21)
# Cyril Matthey-Doret
# 13.10.2016


for bw in 1 5 10;
do
bedtools intersect -a ../data/TAD/merged/short_TADbins"$bw".bed -b ../data/pc_genes/LCL.expressed.pcgene.bed -F 0.25 -wa -c > bin-pcgenes_overlap"$bw"_local.bed;
bedtools intersect -a ../data/TAD/merged/short_TADbins"$bw".bed -b ../data/linc_RNA/LCL.expressed.lincRNA.bed -F 0.25 -wa -c > bin-lincRNA_overlap"$bw"_local.bed;
done

