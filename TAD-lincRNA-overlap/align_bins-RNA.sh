#!bin/bash

# This scripts intend to count the number of overlaps between each TAD bin 
# and at least 25% of the transcripts of protein-coding genes/lincRNAs.
# Cyril Matthey-Doret
# 13.10.2016

module add UHTS/Analysis/BEDTools/2.22.1;

bsub -q priority 'bedtools intersect -a ../data/linc_RNA/LCL.expressed.lincRNA.bed -b ../data/TAD/merged/short_TADbins.bed -f 0.25 -wb > bin-lincRNA_overlap.bed';
bsub -q priority 'bedtools intersect -a ../data/pc_genes/LCL.expressed.pcgenes.bed -b ../data/TAD/merged/short_TADbins.bed -f 0.25 -wb > bin-pcgenes_overlap.bed';

