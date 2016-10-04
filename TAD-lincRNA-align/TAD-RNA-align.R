# This script intends to find lincRNAs with at least 25% of their sequences overlapping the TAD boundaries.
# TAD boundaries are considered to be the regions at the border of TADs, large of 20% the TAD length.


#=====================================================

#Loading data
setwd("/home/cyril/Documents/First_step/data/")
TAD_ori <- read.table("GM12878_TAD_domains.bed")
colnames(TAD_ori) <- c("chr", "start", "end", "gene", "ID")
RNA_ori <- read.table("LCL.expressed.lincRNA.bed")
colnames(RNA_ori) <- c("chr", "start", "end", "gene", "strand")

#====================================================

# creating new dataframes with TAD boundaries
TAD_len <- TAD_ori$end - TAD_ori$start # Computes the length of all RNAs
TAD_boundaries_start <- data.frame(TAD_ori$chr, start= (TAD_ori$start - TAD_len*0.2), end= (TAD_ori$start + TAD_len*0.2), TAD_ori$gene,TAD_ori$ID)
TAD_boundaries_end <- data.frame(TAD_ori$chr, start= (TAD_ori$end - TAD_len*0.2), end= (TAD_ori$end + TAD_len*0.2), TAD_ori$gene, TAD_ori$ID)
# merging dataframes together
TAD_boundaries <- rbind(TAD_boundaries_start,TAD_boundaries_end)

# generating new file with all TAD boundaries.
write.table(TAD_boundaries,file = "TAD_boundaries.bed",sep="\t",quote = F,col.names = F,row.names = F)

# Next step: align lincRNAs to the TAD boundaries. Retain RNAs with at least 25% of their sequence aligned to a TAD boundary.
# calculating length of all RNAs
len_RNA <- cbind(RNA_ori,length=RNA_ori[,3]-RNA_ori[,2])

# Try: aligning -> retrieving alignments -> filtering to keep only alignments of length >25% RNA length.
# see bash script named "alignTAD-RNA.sh" using bedtools intersect.

#====================================================

#Alignment performed, contains duplicates: Some (many) transcripts seem to overlap more than 1 TAD boundary.
overlap_RNA <- read.table("lincRNA_25overlap_TADb.bed")
colnames(overlap_RNA) <- c("chr", "start", "end", "gene", "strand")
length(overlap_RNA$gene[duplicated(overlap_RNA$gene)])
length(overlap_RNA$gene)

#Almost half of the transcript overlap more than 1 TAD boundary. Is this possible ?
#Anyway, removing duplicates
TADbound_lincRNA <- overlap_RNA[!duplicated(overlap_RNA$gene),]
# This means: Out of the 2510 original transcripts, 1420 have at least 25% of their sequence overlapping the TAD boundaries.
write.table(TADbound_lincRNA,file = "TADbound-lincRNA.bed",sep="\t",quote = F,col.names = F,row.names = F)
# Those are the TADbound-lincRNAs
