# This script intends to find lincRNAs and protein coding genes with at least 25% of their sequences overlapping the TAD boundaries.
# TAD boundaries are considered to be the regions at the border of TADs, large of 20% the TAD length.
# Edit: trying different thresholds for definition of TAD boundaries: 5% ,10%, 20%

# Cyril Matthey-Doret, 04.10.2016
options(scipen=999)
#=====================================================

#Loading data
setwd("/home/cyril/Documents/First_step/data/")
TAD_ori <- read.table("TAD/merged/merged_TAD.bed")                #TADs
colnames(TAD_ori) <- c("chr", "start", "end")
RNA_ori <- read.table("linc_RNA/LCL.expressed.lincRNA.bed")              #lincRNAs
colnames(RNA_ori) <- c("chr", "start", "end", "gene", "strand")
pc_ori <- read.table("pc_genes/LCL.expressed.pcgene.bed")                #protein coding genes
colnames(pc_ori) <- c("chr", "start", "end", "gene", "strand")

#====================================================

# creating new dataframes with TAD boundaries.
TAD_len <- TAD_ori$end - TAD_ori$start # Computes the length of all TADs

#5%----------------------------------------------------------

#TADboundaries at the start of every TAD
TAD_boundaries_start <- data.frame(TAD_ori$chr, start= (TAD_ori$start - TAD_len*0.05), end= (TAD_ori$start + TAD_len*0.05))
#Same at the end
TAD_boundaries_end <- data.frame(TAD_ori$chr, start= (TAD_ori$end - TAD_len*0.05), end= (TAD_ori$end + TAD_len*0.05))
# merging dataframes together
TAD_boundaries5 <- rbind(TAD_boundaries_start,TAD_boundaries_end)

#10%--------------------------------------------------------

#TADboundaries at the start of every TAD
TAD_boundaries_start <- data.frame(TAD_ori$chr, start= (TAD_ori$start - TAD_len*0.1), end= (TAD_ori$start + TAD_len*0.1))
#Same at the end
TAD_boundaries_end <- data.frame(TAD_ori$chr, start= (TAD_ori$end - TAD_len*0.1), end= (TAD_ori$end + TAD_len*0.1))
# merging dataframes together
TAD_boundaries10 <- rbind(TAD_boundaries_start,TAD_boundaries_end)

#20%--------------------------------------------------------

#TADboundaries at the start of every TAD
TAD_boundaries_start <- data.frame(TAD_ori$chr, start= (TAD_ori$start - TAD_len*0.2), end= (TAD_ori$start + TAD_len*0.2))
#Same at the end
TAD_boundaries_end <- data.frame(TAD_ori$chr, start= (TAD_ori$end - TAD_len*0.2), end= (TAD_ori$end + TAD_len*0.2))
# merging dataframes together
TAD_boundaries20 <- rbind(TAD_boundaries_start,TAD_boundaries_end)


# generating new file with all TAD boundaries. uncomment to overwrite files.
#write.table(TAD_boundaries5,file = "TAD/merged/TAD_boundaries5.bed",sep="\t",quote = F,col.names = F,row.names = F)
#write.table(TAD_boundaries10,file = "TAD/merged/TAD_boundaries10.bed",sep="\t",quote = F,col.names = F,row.names = F)
#write.table(TAD_boundaries20,file = "TAD/merged/TAD_boundaries20.bed",sep="\t",quote = F,col.names = F,row.names = F)

# Next step: overlap lincRNAs to the TAD boundaries. Retain RNAs with at least 25% of their sequence aligned to a TAD boundary.

# see bash script named "alignTAD-RNA.sh" using bedtools intersect.

#====================================================
#lincRNAs:

#overlaps performed, contains duplicates: Some (many) transcripts seem to overlap more than 1 TAD boundary.
overlap_RNA5 <- read.table("lincRNA_5overlap_TADb.bed")
overlap_RNA10 <- read.table("lincRNA_10overlap_TADb.bed")
overlap_RNA20 <- read.table("lincRNA_20overlap_TADb.bed")
colnames(overlap_RNA5)=colnames(overlap_RNA10)=colnames(overlap_RNA20) <- c("chr", "start", "end", "gene", "strand")

#Number of overlaps: 
length(overlap_RNA5$gene);length(overlap_RNA10$gene);length(overlap_RNA20$gene)

#number of duplicates (lincRNAs matching 2 or more TADb):
length(overlap_RNA5$gene[duplicated(overlap_RNA5$gene)])
length(overlap_RNA10$gene[duplicated(overlap_RNA10$gene)])
length(overlap_RNA20$gene[duplicated(overlap_RNA20$gene)])


#removing duplicates:
TADbound_lincRNA5 <- overlap_RNA5[!duplicated(overlap_RNA5$gene),]
TADbound_lincRNA10 <- overlap_RNA10[!duplicated(overlap_RNA10$gene),]
TADbound_lincRNA20 <- overlap_RNA20[!duplicated(overlap_RNA20$gene),]

#writing into .bed files:
write.table(TADbound_lincRNA5,file = "TADbound-lincRNA5.bed",sep="\t",quote = F,col.names = F,row.names = F)
write.table(TADbound_lincRNA10,file = "TADbound-lincRNA10.bed",sep="\t",quote = F,col.names = F,row.names = F)
write.table(TADbound_lincRNA20,file = "TADbound-lincRNA20.bed",sep="\t",quote = F,col.names = F,row.names = F)

#Identifying nonTADbound-lincRNAs and writing them into files:
nonTADbound_lincRNA5 <- RNA_ori[!(RNA_ori$gene %in% TADbound_lincRNA5$gene),]
write.table(nonTADbound_lincRNA5,file = "nonTADbound-lincRNA5.bed",sep="\t",quote = F,col.names = F,row.names = F)
nonTADbound_lincRNA10 <- RNA_ori[!(RNA_ori$gene %in% TADbound_lincRNA10$gene),]
write.table(nonTADbound_lincRNA10,file = "nonTADbound-lincRNA10.bed",sep="\t",quote = F,col.names = F,row.names = F)
nonTADbound_lincRNA20 <- RNA_ori[!(RNA_ori$gene %in% TADbound_lincRNA20$gene),]
write.table(nonTADbound_lincRNA20,file = "nonTADbound-lincRNA20.bed",sep="\t",quote = F,col.names = F,row.names = F)


#==============================================================
#Protein coding genes:

overlap_pc5 <- read.table("pcgene_5overlap_TADb.bed")
overlap_pc10 <- read.table("pcgene_10overlap_TADb.bed")
overlap_pc20 <- read.table("pcgene_20overlap_TADb.bed")
colnames(overlap_pc5)=colnames(overlap_pc10)=colnames(overlap_pc20) <- c("chr", "start", "end", "gene", "strand")

#number of overlaps found:
length(overlap_pc5$gene); length(overlap_pc10$gene); length(overlap_pc20$gene)

#number of duplicates:
length(overlap_pc5$gene[duplicated(overlap_pc5$gene)])
length(overlap_pc10$gene[duplicated(overlap_pc10$gene)])
length(overlap_pc20$gene[duplicated(overlap_pc20$gene)])

#removing duplicates:
TADbound_pc5 <- overlap_pc5[!duplicated(overlap_pc5$gene),]
TADbound_pc10 <- overlap_pc10[!duplicated(overlap_pc10$gene),]
TADbound_pc20 <- overlap_pc20[!duplicated(overlap_pc20$gene),]

#writing into bed files:
write.table(TADbound_pc5, file = "TADbound-pcgene5.bed", sep="\t", quote=F, col.names = F, row.names = F)
write.table(TADbound_pc10, file = "TADbound-pcgene10.bed", sep="\t", quote=F, col.names = F, row.names = F)
write.table(TADbound_pc20, file = "TADbound-pcgene20.bed", sep="\t", quote=F, col.names = F, row.names = F)
# Those are the TADbound-pcgenes

#Identifying nonTADbound-pcgenes
nonTADbound_pc5 <- pc_ori[!(pc_ori$gene %in% TADbound_pc5$gene),]
write.table(nonTADbound_pc5,file = "nonTADbound-pcgene5.bed",sep="\t",quote = F,col.names = F,row.names = F)

nonTADbound_pc10 <- pc_ori[!(pc_ori$gene %in% TADbound_pc10$gene),]
write.table(nonTADbound_pc10,file = "nonTADbound-pcgene10.bed",sep="\t",quote = F,col.names = F,row.names = F)

nonTADbound_pc20 <- pc_ori[!(pc_ori$gene %in% TADbound_pc20$gene),]
write.table(nonTADbound_pc20,file = "nonTADbound-pcgene20.bed",sep="\t",quote = F,col.names = F,row.names = F)


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

#==================================================
# Dividing TAD boundaries into bins 
setwd("/Users/cmatthe5/Documents/First_step/data")
nTADb_lincRNA10 <- read.table("linc_RNA/nonTADbound-lincRNA10.bed")
TADb_lincRNA10 <- read.table("linc_RNA/nonTADbound-lincRNA10.bed")
TADb_pcgene10 <- read.table("pc_genes/TADbound-pcgene10.bed")
nTADb_pcgene10 <- read.table("pc_genes/nonTADbound-pcgene10.bed")
colnames(nTADb_lincRNA10)=colnames(TADb_pcgene10)=
  colnames(TADb_lincRNA10) = colnames(nTADb_pcgene10)<-c("chr", "start","end","ID", "strand")
TADb10 <- read.table("TAD/TAD_boundaries10.bed")
TADb10 <- TADb10[,-4]
TAD <- read.table("TAD/GM12878_TAD_domains.bed")
TAD <- TAD[,-4]
colnames(TADb10) = colnames(TAD)<-c("chr", "start","end","ID")

TADbins <- 
