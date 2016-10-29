# This script tests for enrichment in enhancer-bound and non-enhancer bound protein-coding genes and lincRNAs in TAD-bins.
# This is done by a series of GAT calls using different values for segment, workspace and annotation arguments.
# Cyril Matthey-Doret
# 29.10.2016
g = "genes/"
b = "bins/"
count = 1
W = "hg19.genome.bed.gat" 
A = '$b''short_bins5.bed'
for S in '$g''elinc_pr.bed' '$g''elinc_prb.bed' '$g''epc_pr.bed' '$g''epc_prb.bed' '$g''nelinc_pr.bed' '$g''nelinc_prb.bed' '$g''nepc_pr.bed' '$g''nepc_prb.bed'
do
    gat-run.py  --verbose=5\
                --log="$count".log\
                --segments=S\
                --annotation=A\
                --workspace=W\
    count = $((count+1))           
    gat-run.py  --verbose=5\
                --log="r""$count".log\
                --segments=A\
                --annotation=S\
                --workspace=W\
    count = $((count+1))           
done


W = "intergenic.bed"
count = 1
for S in 
do
    for A in 
    do

    gat-run.py  --verbose=5\
                --log="$count".log\
                --segments=S\
                --annotation=A\
                --workspace=W\
    count = $((count+1))           
    done
done
