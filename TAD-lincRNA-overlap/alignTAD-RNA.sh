#!/bin/bash

#This script intends to align lincRNA and protein coding genes with TAD boundaries. It outputs RNAs which have at least
#25% of their sequence overlapping a TAD boundary.
# Cyril Matthey-Doret, 04.10.2016

bsub -q priority 'bedtools intersect -a ../data/LCL.expressed.lincRNA.bed -b ../data/TAD_boundaries5.bed -f 0.25 -wa > lincRNA_5overlap_TADb.bed';

bsub -q priority 'bedtools intersect -a ../data/LCL.expressed.pcgene.bed -b ../data/TAD_boundaries5.bed -f 0.25 -wa > pcgene_5overlap_TADb.bed';

bsub -q priority 'bedtools intersect -a ../data/LCL.expressed.lincRNA.bed -b ../data/TAD_boundaries10.bed -f 0.25 -wa > lincRNA_10overlap_TADb.bed';

bsub -q priority 'bedtools intersect -a ../data/LCL.expressed.pcgene.bed -b ../data/TAD_boundaries10.bed -f 0.25 -wa > pcgene_10overlap_TADb.bed';

bsub -q priority 'bedtools intersect -a ../data/LCL.expressed.lincRNA.bed -b ../data/TAD_boundaries25.bed -f 0.25 -wa > lincRNA_20overlap_TADb.bed';

bsub -q priority 'bedtools intersect -a ../data/LCL.expressed.pcgene.bed -b ../data/TAD_boundaries25.bed -f 0.25 -wa > pcgene_20overlap_TADb.bed'
