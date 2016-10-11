#!/bin/bash

#This script intends to align lincRNA and protein coding genes with TAD boundaries of merged TADs. It outputs RNAs which have at least
#25% of their sequence overlapping a TAD boundary.
# Cyril Matthey-Doret, 04.10.2016

for t in 5 10 20;
do
bsub -q priority 'bedtools intersect -a ../data/LCL.expressed.lincRNA.bed -b ../data/TAD_boundaries${t}.bed -f 0.25 -wa > lincRNA_${t}overlap_TADb.bed';

bsub -q priority 'bedtools intersect -a ../data/LCL.expressed.pcgene.bed -b ../data/TAD_boundaries${t}.bed -f 0.25 -wa > pcgene_${t}overlap_TADb.bed';

done
