# This script intends to find lincRNAs and protein coding genes with at least 25% of their sequences overlapping the TAD boundaries.
# TAD boundaries are considered to be the regions at the border of TADs, large of 20% the TAD length.

# Cyril Matthey-Doret, 04.10.2016

#=====================================================

#Loading data
setwd("/home/cyril/Documents/First_step/data/")
TAD_ori <- read.table("GM12878_TAD_domains.bed")                #TADs
colnames(TAD_ori) <- c("chr", "start", "end", "gene", "ID")
RNA_ori <- read.table("LCL.expressed.lincRNA.bed")              #lincRNAs
colnames(RNA_ori) <- c("chr", "start", "end", "gene", "strand")
pc_ori <- read.table("LCL.expressed.pcgene.bed")                #protein coding genes
colnames(pc_ori) <- c("chr", "start", "end", "gene", "strand")

#====================================================

# creating new dataframes with TAD boundaries
TAD_len <- TAD_ori$end - TAD_ori$start # Computes the length of all RNAs
#TADboundaries at the start of every TAD
TAD_boundaries_start <- data.frame(TAD_ori$chr, start= (TAD_ori$start - TAD_len*0.2), end= (TAD_ori$start + TAD_len*0.2), TAD_ori$gene,TAD_ori$ID)
#Same at the end
TAD_boundaries_end <- data.frame(TAD_ori$chr, start= (TAD_ori$end - TAD_len*0.2), end= (TAD_ori$end + TAD_len*0.2), TAD_ori$gene, TAD_ori$ID)
# merging dataframes together
TAD_boundaries <- rbind(TAD_boundaries_start,TAD_boundaries_end)

# generating new file with all TAD boundaries.
write.table(TAD_boundaries,file = "TAD_boundaries.bed",sep="\t",quote = F,col.names = F,row.names = F)

# Next step: overlap lincRNAs to the TAD boundaries. Retain RNAs with at least 25% of their sequence aligned to a TAD boundary.

# see bash script named "alignTAD-RNA.sh" using bedtools intersect.

#====================================================
#lincRNAs:

#overlap performed, contains duplicates: Some (many) transcripts seem to overlap more than 1 TAD boundary.
overlap_RNA <- read.table("lincRNA_25overlap_TADb.bed")
colnames(overlap_RNA) <- c("chr", "start", "end", "gene", "strand")


length(overlap_RNA$gene)
# lincRNA: 2554 overlap events were found by the intersect program with -f 0.25

length(overlap_RNA$gene[duplicated(overlap_RNA$gene)])
# lincRNA: 1134 are duplicates

#Anyway, removing duplicates
TADbound_lincRNA <- overlap_RNA[!duplicated(overlap_RNA$gene),]
# lincRNA: Out of the 2510 original transcripts, 1420 have at least 25% of their sequence overlapping the TAD boundaries.
write.table(TADbound_lincRNA,file = "TADbound-lincRNA.bed",sep="\t",quote = F,col.names = F,row.names = F)
# Those are the TADbound-lincRNAs

#Identifying nonTADbound-lincRNAs
nonTADbound_lincRNA <- RNA_ori[!(RNA_ori$gene %in% TADbound_lincRNA$gene),]
write.table(nonTADbound_lincRNA,file = "nonTADbound-lincRNA.bed",sep="\t",quote = F,col.names = F,row.names = F)
#Those are the nonTADbound-lincRNAs (there are 1090 of these)

#==============================================================
#Protein coding genes:

overlap_pc <- read.table("pcgene_25overlap_TADb.bed")
colnames(overlap_pc) <- c("chr", "start", "end", "gene", "strand")

length(overlap_pc$gene)
# pcgenes: 21544 overlap events were found using the same parameters.

length(overlap_pc$gene[duplicated(overlap_pc$gene)])
# pcgenes: 10207 are duplicates

TADbound_pc <- overlap_pc[!duplicated(overlap_pc$gene),]
#pcgenes: Out of the 14846 original genes, 11337 have at least 25% of their sequence overlapping the TAD boundaries.
write.table(TADbound_pc, file = "TADbound-pcgene.bed", sep="\t", quote=F, col.names = F, row.names = F)
# Those are the TADbound-pcgenes

#Identifying nonTADbound-pcgenes
nonTADbound_pc <- pc_ori[!(pc_ori$gene %in% TADbound_pc$gene),]
write.table(nonTADbound_pc,file = "nonTADbound-pcgene.bed",sep="\t",quote = F,col.names = F,row.names = F)
#Those are the nonTADbound-pcgenes (there are 3509 of these)

#===================================================
# Visualizing length of lincRNAs, pcgenes and TAD boundaries.
par(mfrow=c(1,3))
hist(log10(TAD_boundaries$end - TAD_boundaries$start))
hist(log10(TADbound_lincRNA$end - TADbound_lincRNA$start))
hist(log10(RNA_ori$end - RNA_ori$start))
hist(log10(TADbound_pc$end - TADbound_pc$start))
hist(log10(pc_ori$end - pc_ori$start))
#Summary statistics
summary(TAD_boundaries$end - TAD_boundaries$start)
#lincRNA:
summary(TADbound_lincRNA$end - TADbound_lincRNA$start)
summary(RNA_ori$end - RNA_ori$start)
summary(nonTADbound_lincRNA$end - nonTADbound_lincRNA$start)
#pcgenes
summary(TADbound_pc$end - TADbound_pc$start)
summary(pc_ori$end - pc_ori$start)
summary(nonTADbound_pc$end - nonTADbound_pc$start)


