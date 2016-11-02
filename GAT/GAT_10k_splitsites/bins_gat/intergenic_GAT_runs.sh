#!/bin/bash

#BSUB -J inter_gat
#BSUB -M 16777216
# This script tests for enrichment in enhancer-bound and non-enhancer bound lincRNAs in TAD-bins.
# This is done by a series of GAT calls using different values for segment and annotation arguments in intergenic space of the genome.
# Cyril Matthey-Doret
# 29.10.2016


g=$(find ../data/GAT/genes/*linc*)
b="../data/GAT/bins/short_bins10.bed"
W="../data/GAT/intergenic.bed" 
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
                --num-samples=10000 \
                --qvalue-method=BH \
                --ignore-segment-tracks \
                --isochore-file="../data/GAT/hg19.fa.corr_term_ISOisochore.bed" \
                >'gat_'$desc'.tsv'
    
    desc="W_intergenic_S_short10bins_A_""${sS%\.*}"
    gat-run.py  --verbose=5 \
                --log='log_'$desc'.log' \
                --num-samples=10000 \
                --qvalue-method=BH \
                --segment-file=$A \
                --annotation-file=$S \
                --workspace=$W \
                --isochore-file="../data/GAT/hg19.fa.corr_term_ISOisochore.bed" \
                >'gat_'$desc'.tsv'        
done


 
