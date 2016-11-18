# In this script I compare the sequence conservation of enhancer-associated and promoter-associated lincRNAs and protein-coding genes. 
# This is done by comparing the distribution of their average phastCons scores (1 score per nucleotide, averaged over each gene).
# I also use the distribution of ancestral repeats phastCons scores as a comparison, since those repeats are assumed to be evolving
# neutrally and were already present before the divergence of human and mouse. Genes that have a higher phastCons scores compared to ARs
# are likely to be under purifying selection, whereas those with lower phastCons scores should be under positive selection

# Cyril Matthey-Doret
# 25.10.2016
##################################

# Loading data:

setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
setwd("/Users/cmatthe5/Documents/First_step/data/")

elinc <- read.table("enhancer_bound/elinc_prb.bed")
plinc <- read.table("enhancer_bound/plinc_prb.bed")
epc <- read.table("enhancer_bound/epc_prb.bed")
ppc <- read.table("enhancer_bound/ppc_prb.bed")
pri.cons_linc <- read.table("seq_conserv/all.lincRNA.primates.phascons.exons.txt", header = T)
pri.cons_pc <- read.table("seq_conserv/all.pcgene.primates.phascons.exons.txt", header = T)
pri.cons_AR <- read.table("seq_conserv/AR.primates.phascons.txt", header = T)
mam.cons_linc <- read.table("seq_conserv/all.lincRNA.mammals.phascons.exons.txt",header = T)
mam.cons_pc <- read.table("seq_conserv/all.pcgene.mammals.phascons.exons.txt",header = T)
mam.cons_AR <- read.table("seq_conserv/AR.mammals.phascons.txt", header=T)

colnames(elinc)=colnames(plinc)=colnames(epc)=colnames(ppc) <-c("chr", "start", "end", "gene", "strand")

#================================

# Splitting prom.assoc vs enh.assoc

mam.cons_elinc <-mam.cons_linc[mam.cons_linc$gene %in% elinc$gene,]
mam.cons_plinc <-mam.cons_linc[mam.cons_linc$gene %in% plinc$gene,]
mam.cons_epc <-mam.cons_pc[mam.cons_pc$gene %in% epc$gene,]
mam.cons_ppc <-mam.cons_pc[mam.cons_pc$gene %in% ppc$gene,]

pri.cons_elinc <-pri.cons_linc[pri.cons_linc$gene %in% elinc$gene,]
pri.cons_plinc <-pri.cons_linc[pri.cons_linc$gene %in% plinc$gene,]
pri.cons_epc <-pri.cons_pc[pri.cons_pc$gene %in% epc$gene,]
pri.cons_ppc <-pri.cons_pc[pri.cons_pc$gene %in% ppc$gene,]

whole_cons <- rbind(cbind(mam.cons_elinc, gr=rep("mam"),assoc=rep("e"),gentype=rep("lincRNA")),
      cbind(mam.cons_plinc,gr=rep("mam"),assoc=rep("p"),gentype=rep("lincRNA")),
      cbind(mam.cons_epc,gr=rep("mam"),assoc=rep("e"),gentype=rep("pc")),
      cbind(mam.cons_ppc,gr=rep("mam"),assoc=rep("p"),gentype=rep("pc")),
      cbind(pri.cons_elinc,gr=rep("pri"),assoc=rep("e"),gentype=rep("lincRNA")),
      cbind(pri.cons_plinc,gr=rep("pri"),assoc=rep("p"),gentype=rep("lincRNA")),
      cbind(pri.cons_epc,gr=rep("pri"),assoc=rep("e"),gentype=rep("pc")),
      cbind(pri.cons_ppc,gr=rep("pri"),assoc=rep("p"),gentype=rep("pc")),
      cbind(mam.cons_AR,gr=rep("mam"),assoc=rep("AR"),gentype=rep("AR")),
      cbind(pri.cons_AR,gr=rep("pri"),assoc=rep("AR"),gentype=rep("AR")))
write.table(x = whole_cons,file = "seq_conserv/enhancer_bound/whole_cons.txt",quote = F,sep = "\t",row.names = F,col.names = T)

# Visualisations are directly implemented in report 5.
