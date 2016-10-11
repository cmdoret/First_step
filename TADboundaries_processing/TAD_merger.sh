# This script intends to merge overlapping TADs together, in order to consider only the large ones.
# This should reduce the noise caused by the overlappiing boundaries of smaller TADs nested inside the larger ones.
# Cyril Matthey-Doret
# 11.10.2016


bsub -q priority 'bedtools merge -i ../data/TAD/GM12878_TAD_domains.bed > merged_TAD.bed'

