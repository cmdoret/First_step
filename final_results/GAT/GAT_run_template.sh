#!/bin/bash

# This script is a template for testing for enrichment of a desired set of genes, peaks or any other segment in a desired set of annotations.
# This is done by a series of GAT calls using different values for segment  arguments in given workspace.
# Cyril Matthey-Doret
# 29.10.2016

g=$(find ../path/*pattern*)
b="../path/annotation"
W="../path/workspace" 
A=$b
for S in $g # Iterating over segments files matching desired pattern in the folder
do
    sS=${S##*/} # removing path from the segment filename
    desc="W_<workspace>_S_"${sS%.*}"_A_<annotation>" # Removing extension from segment filename and giving meaningful filename to the future output.
    gat-run.py  --verbose=5 \
                --log='log_'$desc'.log' \
                --segment-file=$S \
                --annotation-file=$A \
                --workspace=$W \
                --ignore-segment-tracks \
                --num-samples=10000 \
                --counter=segment-overlap \
                --qvalue-method=BH \
                --isochore-file="../../data/GAT/hg19.fa.corr_term_ISOisochore.bed" \
                >'gat_'$desc'.tsv'
done    

