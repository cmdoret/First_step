# In this script I compare the sequence conservation of TADbound and non-TADbound lincRNAs and protein-coding genes. 
# This is done by comparing the distribution of their average phastCons scores (1 score per nucleotide, averaged over each gene).
# I also use the distribution of ancestral repeats phastCons scores as a comparison, since those repeats are assumed to be evolving
# neutrally and were already present before the divergence of human and mouse. Genes that have a higher phastCons scores compared to ARs
# are likely to be under purifying selection, whereas those with lower phastCons scores should be under positive selection

# Cyril Matthey-Doret
# 25.10.2016
##################################

# Loading data:

setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")

TADb_linc <- read.table("linc_RNA/merged/flexTADbound-lincRNA5.bed")
nTADb_linc <- read.table("linc_RNA/merged/flexnonTADbound-lincRNA5.bed")
TADb_pc <- read.table("pc_genes/merged/flexTADbound-pcgene5.bed")
nTADb_pc <- read.table("pc_genes/merged/flexnonTADbound-pcgene5.bed")
pri.cons_linc <- read.table("seq_conserv/all.lincRNA.primates.phascons.exons.txt", header = T)
pri.cons_pc <- read.table("seq_conserv/all.pcgene.primates.phascons.exons.txt", header = T)
pri.cons_AR <- read.table("seq_conserv/AR.primates.phascons.txt", header = T)
mam.cons_linc <- read.table("seq_conserv/all.lincRNA.mammals.phascons.exons.txt",header = T)
mam.cons_pc <- read.table("seq_conserv/all.pcgene.mammals.phascons.exons.txt",header = T)
mam.cons_AR <- read.table("seq_conserv/AR.mammals.phascons.txt", header=T)

colnames(TADb_linc)=colnames(TADb_pc)=colnames(nTADb_linc)=colnames(nTADb_pc) <-c("chr", "start", "end", "gene", "strand")

#================================

# Splitting nTADb vs TADb

mam.cons_Tb_linc <-mam.cons_linc[mam.cons_linc$gene %in% TADb_linc$gene,]
mam.cons_nTb_linc <-mam.cons_linc[mam.cons_linc$gene %in% nTADb_linc$gene,]
mam.cons_Tb_pc <-mam.cons_pc[mam.cons_pc$gene %in% TADb_pc$gene,]
mam.cons_nTb_pc <-mam.cons_linc[mam.cons_linc$gene %in% nTADb_linc$gene,]

pri.cons_Tb_linc <-pri.cons_linc[pri.cons_linc$gene %in% TADb_linc$gene,]
pri.cons_nTb_linc <-pri.cons_linc[pri.cons_linc$gene %in% nTADb_linc$gene,]
pri.cons_Tb_pc <-pri.cons_pc[pri.cons_pc$gene %in% TADb_pc$gene,]
pri.cons_nTb_pc <-pri.cons_pc[pri.cons_pc$gene %in% nTADb_pc$gene,]

whole_cons <- rbind(cbind(mam.cons_Tb_linc, gr=rep("mam"),TAD=rep("TADb"),gentype=rep("linc")),
      cbind(mam.cons_nTb_linc,gr=rep("mam"),TAD=rep("nTADb"),gentype=rep("linc")),
      cbind(mam.cons_Tb_pc,gr=rep("mam"),TAD=rep("TADb"),gentype=rep("pc")),
      cbind(mam.cons_nTb_pc,gr=rep("mam"),TAD=rep("nTADb"),gentype=rep("pc")),
      cbind(pri.cons_Tb_linc,gr=rep("pri"),TAD=rep("TADb"),gentype=rep("linc")),
      cbind(pri.cons_nTb_linc,gr=rep("pri"),TAD=rep("nTADb"),gentype=rep("linc")),
      cbind(pri.cons_Tb_pc,gr=rep("pri"),TAD=rep("TADb"),gentype=rep("pc")),
      cbind(pri.cons_nTb_pc,gr=rep("pri"),TAD=rep("nTADb"),gentype=rep("pc")),
      cbind(mam.cons_AR,gr=rep("mam"),TAD=rep("AR"),gentype=rep("AR")),
      cbind(pri.cons_AR,gr=rep("pri"),TAD=rep("AR"),gentype=rep("AR")))
write.table(x = whole_cons,file = "seq_conserv/flex_whole_cons.txt",quote = F,sep = "\t",row.names = F,col.names = T)

# Visualisations are directly implemented in report 5.
