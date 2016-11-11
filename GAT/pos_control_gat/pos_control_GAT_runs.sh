#!/bin/bash

#BSUB -M 16777216
#BSUB -J whole_gat
# This script tests for enrichment of transcription factors binding sites in TAD bins
# This is odone using GAT with whole genome as the workspace,binding sites as the segments
# and TAD bins as annotations
# Cyril Matthey-Doret
# 29.10.2016


g="../data/chip_seq/GM12878*"
b="../data/GAT/bins/short_bins10.bed"
W="../data/GAT/hg19.genome.bed.gat" 
A=$b
for S in $g
do
    sS=${S##*/}
	    desc="W_intergenic_S_"${sS%.*}"_A_short10bins"
	    gat-run.py  --verbose=5 \
	                --log='log_'$desc'.log' \
	                --segment-file=$S \
	                --annotation-file=$A \
	                --workspace=$W \
	                --counter=segment-overlap \
	                --ignore-segment-tracks \
	                --num-samples=10000 \
	                --qvalue-method=BH \
	                --isochore-file="../../data/GAT/hg19.fa.corr_term_ISOisochore.bed" \
	                >'gat_'$desc'.tsv'
done


