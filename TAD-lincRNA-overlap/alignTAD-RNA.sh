#!/bin/bash

#This script intends to align lincRNA and protein coding genes with TAD boundaries. It outputs RNAs which have at least
#25% of their sequence overlapping a TAD boundary.

bsub -q priority 'bedtools intersect -a ../data/LCL.expressed.lincRNA.bed -b ../data/GM12878_TAD_boundaries.bed -f 0.25 -wa > lincRNA_25overlap_TADb.bed';

bsub -q priority 'bedtools intersect -a ../data/LCL.expressed.pcgene.bed -b ../data/GM12878_TAD_boundaries.bed -f 0.25 -wa > pcgene_25overlap_TADb.bed'
