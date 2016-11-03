#!/bin/bash
# This script tests for enrichment in enhancer-bound and non-enhancer bound protein-coding genes and lincRNAs in TAD-bins.
# This is done by a series of GAT calls using different values for segment and annotation arguments in the whole genome.
# Cyril Matthey-Doret
# 29.10.2016


g=$(find ../../data/GAT/genes/elinc*)
b="../../data/GAT/bins/short_bins10.bed"
W="../../data/linc_RNA/expressed_linc_workspace.bed" 
A=$b
for S in $g
do
    sS=${S##*/}
    desc="W_lincRNA_S_"${sS%.*}"_A_short10bins"
    gat-run.py  --verbose=5 \
                --log='log_'$desc'.log' \
                --segment-file=$S \
                --annotation-file=$A \
                --workspace=$W \
                --ignore-segment-tracks \
                --num-samples=10000 \
                --qvalue-method=BH \
                --isochore-file="../../data/GAT/hg19.fa.corr_term_ISOisochore.bed" \
                >'gat_'$desc'.tsv'
    
    desc="W_lincRNA_S_short10bins_A_""${sS%\.*}"
    gat-run.py  --verbose=5 \
                --log='log_'$desc'.log' \
                --num-samples=10000 \
                --qvalue-method=BH \
                --segment-file=$A \
                --counter=segment-overlap \
                --annotation-file=$S \
                --workspace=$W \
                --isochore-file="../../data/GAT/hg19.fa.corr_term_ISOisochore.bed" \
                >'gat_'$desc'.tsv'        
done


