# In this script I check if there are correlations between the expression of TAD-bound lincRNAs and different properties of these lincRNAs.
# Cyril Matthey-Doret
# Tue Oct 25 22:32:48 2016 ------------------------------

#Loading data:
#setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
setwd("/home/cyril/Documents/First_step/data/")
TAD <- read.table("TAD/merged/merged_TAD.bed")
colnames(TAD)<-c("chr","start","end")
Tb_linc <- read.table("linc_RNA/merged/flexTADbound-lincRNA5.bed")
Tb_pc <- read.table("pc_genes/merged/flexTADbound-pcgene5.bed")
colnames(Tb_pc)=colnames(Tb_linc) <- c("chr","start","end","gene","strand","TAD_chr","TAD_start","TAD_end","TAD_ID")
nTb_linc <- read.table("linc_RNA/merged/flexnonTADbound-lincRNA5.bed")
nTb_pc <- read.table("pc_genes/merged/flexnonTADbound-pcgene5.bed")
colnames(nTb_pc)=colnames(nTb_linc) <- c("chr","start","end","gene","strand")
TAD$ID <- paste0("merged_",seq(1:length(TAD[,1])))
exp_Tb_linc <- read.table("expression/merged/flex_exp_Tb_lincRNA5.txt")
exp_nTb_linc <-read.table("expression/merged/flex_exp_nTb_lincRNA5.txt") 
exp_Tb_pc <- read.table("expression/merged/flex_exp_Tb_pc5.txt")
exp_nTb_pc <-read.table("expression/merged/flex_exp_nTb_pc5.txt") 
colnames(exp_Tb_linc)=colnames(exp_nTb_linc)=colnames(exp_Tb_pc)=colnames(exp_nTb_pc) <- c("gene","expression")

#=======================================================
# Expression vs length of TAD:

Tb_linc <- merge(Tb_linc, exp_Tb_linc,by="gene", all=F)
nTb_linc <- merge(nTb_linc, exp_nTb_linc,by="gene", all=F)
Tb_pc <- merge(Tb_pc, exp_Tb_pc,by="gene", all=F)
nTb_pc <- merge(nTb_pc, exp_nTb_pc,by="gene", all=F)

smoothScatter(log10(Tb_linc$expression),log10(Tb_linc$TAD_end - Tb_linc$TAD_start))
smoothScatter(log10(Tb_pc$expression),log10(Tb_pc$TAD_end - Tb_pc$TAD_start))
cor.test(Tb_linc$expression,(Tb_linc$TAD_end - Tb_linc$TAD_start))
cor.test(Tb_pc$expression,(Tb_pc$TAD_end - Tb_pc$TAD_start))
# No correlation

#======================================================
# Expression vs number of pcgenes in the same TAD:

read.table("")
TAD$countpc <- rep(0)
for(t in row(TAD)){
  TAD$countpc[t] <- length(Tb_pc$gene[Tb_pc$TAD_ID==TAD$ID[t]])
}
Tb_linc$countpc <- rep(0)
for(l in row(Tb_linc)){
  Tb_linc$countpc[l] <-TAD$countpc[TAD$ID==Tb_linc$TAD_ID[l]]
}
smoothScatter(Tb_linc$expression,)
#Question: What is "in a TAD"? in:
  #a) same TAD boundary
  #b) same TAD boundary or other boundary of same TAD
  #c) same TAD sensus stricto
  #d) same TAD + TAD boundaries.

# Expression vs median expression levels of pcgenes in the same TAD:

# Expression vs amount of chromosomal interaction (Hi-C data):

