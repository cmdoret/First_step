# The purpose of this script is to generate a summarized table containing PCG and lincRNAs with their associated sequence conservation
# (average phastCons scores) throughout both mammalian and primate evolution, and their category based on promoter/enhancer overlap.

# Cyril Matthey-Doret
# 25.10.2016
##################################

# Loading data:
setwd("/Users/cmatthe5/Documents/First_step/data/")
setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")

# Loading lincRNAs sets
ne_linc <- read.table("enhancer_bound/all_combinations/ne_linc_pr.bed")  # Overlap no enhancer marks, does not take promoters into account
e_linc <- read.table("enhancer_bound/all_combinations/e_linc_pr.bed")  # Overlap enhancer marks, does not take promoters into account
ne.p_linc <- read.table("enhancer_bound/all_combinations/ne.p_linc_pr.bed")  # Overlap promoter marks, but no enhancer marks
e.np_linc <- read.table("enhancer_bound/all_combinations/e.np_linc_pr.bed")  # Overlap enhancer marks, but no promoter marks
e.p_linc <- read.table("enhancer_bound/all_combinations/e.p_linc_pr.bed")  # Overlap promoter marks, and enhancer marks
ne.np_linc <- read.table("enhancer_bound/all_combinations/ne.np_linc_pr.bed")  # Overlap neither promoter marks, nor enhancer marks

# Loading pcgenes sets
ne_pc <- read.table("enhancer_bound/all_combinations/ne_pc_pr.bed")  # Overlap no enhancer marks, does not take promoters into account
e_pc <- read.table("enhancer_bound/all_combinations/e_pc_pr.bed")  # Overlap enhancer marks, does not take promoters into account
ne.p_pc <- read.table("enhancer_bound/all_combinations/ne.p_pc_pr.bed")  # Overlap promoter marks, but no enhancer marks
e.np_pc <- read.table("enhancer_bound/all_combinations/e.np_pc_pr.bed")  # Overlap enhancer marks, but no promoter marks
e.p_pc <- read.table("enhancer_bound/all_combinations/e.p_pc_pr.bed")  # Overlap promoter marks, and enhancer marks
ne.np_pc <- read.table("enhancer_bound/all_combinations/ne.np_pc_pr.bed")  # Overlap neither promoter marks, nor enhancer marks

pri.cons_lincRNA <- read.table("seq_conserv/all.lincRNA.primates.phascons.exons.txt", header = T)
pri.cons_pcgene <- read.table("seq_conserv/all.pcgene.primates.phascons.exons.txt", header = T)
pri.cons_AR <- read.table("seq_conserv/AR.primates.phascons.txt", header = T)
mam.cons_lincRNA <- read.table("seq_conserv/all.lincRNA.mammals.phascons.exons.txt",header = T)
mam.cons_pcgene <- read.table("seq_conserv/all.pcgene.mammals.phascons.exons.txt",header = T)
mam.cons_AR <- read.table("seq_conserv/AR.mammals.phascons.txt", header=T)

colnames(ne_linc)= colnames(e_linc)=colnames(ne.p_linc)=colnames(e.np_linc)=colnames(e.p_linc)=
  colnames(ne.np_linc)=colnames(ne_pc)=colnames(e_pc)=colnames(ne.p_pc)=colnames(e.np_pc)=colnames(e.p_pc)=
  colnames(ne.np_pc) <-c("chr", "start", "end", "gene", "strand")

#================================
#Splitting conservation data into categories
mam.cons_ne_linc <-mam.cons_lincRNA[mam.cons_lincRNA$gene %in% ne_linc$gene,]
mam.cons_e_linc <-mam.cons_lincRNA[mam.cons_lincRNA$gene %in% e_linc$gene,]
mam.cons_e.p_linc <-mam.cons_lincRNA[mam.cons_lincRNA$gene %in% e.p_linc$gene,]
mam.cons_ne.p_linc <-mam.cons_lincRNA[mam.cons_lincRNA$gene %in% ne.p_linc$gene,]
mam.cons_e.np_linc <-mam.cons_lincRNA[mam.cons_lincRNA$gene %in% e.np_linc$gene,]
mam.cons_ne.np_linc <-mam.cons_lincRNA[mam.cons_lincRNA$gene %in% ne.np_linc$gene,]

mam.cons_ne_pc <-mam.cons_pcgene[mam.cons_pcgene$gene %in% ne_pc$gene,]
mam.cons_e_pc <-mam.cons_pcgene[mam.cons_pcgene$gene %in% e_pc$gene,]
mam.cons_e.p_pc <-mam.cons_pcgene[mam.cons_pcgene$gene %in% e.p_pc$gene,]
mam.cons_ne.p_pc <-mam.cons_pcgene[mam.cons_pcgene$gene %in% ne.p_pc$gene,]
mam.cons_e.np_pc <-mam.cons_pcgene[mam.cons_pcgene$gene %in% e.np_pc$gene,]
mam.cons_ne.np_pc <-mam.cons_pcgene[mam.cons_pcgene$gene %in% ne.np_pc$gene,]

pri.cons_ne_linc <-pri.cons_lincRNA[pri.cons_lincRNA$gene %in% ne_linc$gene,]
pri.cons_e_linc <-pri.cons_lincRNA[pri.cons_lincRNA$gene %in% e_linc$gene,]
pri.cons_e.p_linc <-pri.cons_lincRNA[pri.cons_lincRNA$gene %in% e.p_linc$gene,]
pri.cons_ne.p_linc <-pri.cons_lincRNA[pri.cons_lincRNA$gene %in% ne.p_linc$gene,]
pri.cons_e.np_linc <-pri.cons_lincRNA[pri.cons_lincRNA$gene %in% e.np_linc$gene,]
pri.cons_ne.np_linc <-pri.cons_lincRNA[pri.cons_lincRNA$gene %in% ne.np_linc$gene,]

pri.cons_ne_pc <-pri.cons_pcgene[pri.cons_pcgene$gene %in% ne_pc$gene,]
pri.cons_e_pc <-pri.cons_pcgene[pri.cons_pcgene$gene %in% e_pc$gene,]
pri.cons_e.p_pc <-pri.cons_pcgene[pri.cons_pcgene$gene %in% e.p_pc$gene,]
pri.cons_ne.p_pc <-pri.cons_pcgene[pri.cons_pcgene$gene %in% ne.p_pc$gene,]
pri.cons_e.np_pc <-pri.cons_pcgene[pri.cons_pcgene$gene %in% e.np_pc$gene,]
pri.cons_ne.np_pc <-pri.cons_pcgene[pri.cons_pcgene$gene %in% ne.np_pc$gene,]

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

whole_cons <- rbind(cbind(mam.cons_e_linc,promoter=rep("both"),enhancer=rep("+"),gentype=rep("linc"),gr=rep("mam")), 
                   cbind(mam.cons_ne_linc,promoter=rep("both"),enhancer=rep("-"),gentype=rep("linc"),gr=rep("mam")),
                   cbind(mam.cons_e.np_linc,promoter=rep("-"),enhancer=rep("+"),gentype=rep("linc"),gr=rep("mam")),
                   cbind(mam.cons_ne.p_linc,promoter=rep("+"),enhancer=rep("-"),gentype=rep("linc"),gr=rep("mam")),
                   cbind(mam.cons_e.p_linc,promoter=rep("+"),enhancer=rep("+"),gentype=rep("linc"),gr=rep("mam")),
                   cbind(mam.cons_ne.np_linc,promoter=rep("-"),enhancer=rep("-"),gentype=rep("linc"),gr=rep("mam")),
                   cbind(mam.cons_e_pc,promoter=rep("both"),enhancer=rep("+"),gentype=rep("pc"),gr=rep("mam")), 
                   cbind(mam.cons_ne_pc,promoter=rep("both"),enhancer=rep("-"),gentype=rep("pc"),gr=rep("mam")),
                   cbind(mam.cons_e.np_pc,promoter=rep("-"),enhancer=rep("+"),gentype=rep("pc"),gr=rep("mam")),
                   cbind(mam.cons_ne.p_pc,promoter=rep("+"),enhancer=rep("-"),gentype=rep("pc"),gr=rep("mam")),
                   cbind(mam.cons_e.p_pc,promoter=rep("+"),enhancer=rep("+"),gentype=rep("pc"),gr=rep("mam")),
                   cbind(mam.cons_ne.np_pc,promoter=rep("-"),enhancer=rep("-"),gentype=rep("pc"),gr=rep("mam")),
                   cbind(pri.cons_e_linc,promoter=rep("both"),enhancer=rep("+"),gentype=rep("linc"),gr=rep("pri")), 
                   cbind(pri.cons_ne_linc,promoter=rep("both"),enhancer=rep("-"),gentype=rep("linc"),gr=rep("pri")),
                   cbind(pri.cons_e.np_linc,promoter=rep("-"),enhancer=rep("+"),gentype=rep("linc"),gr=rep("pri")),
                   cbind(pri.cons_ne.p_linc,promoter=rep("+"),enhancer=rep("-"),gentype=rep("linc"),gr=rep("pri")),
                   cbind(pri.cons_e.p_linc,promoter=rep("+"),enhancer=rep("+"),gentype=rep("linc"),gr=rep("pri")),
                   cbind(pri.cons_ne.np_linc,promoter=rep("-"),enhancer=rep("-"),gentype=rep("linc"),gr=rep("pri")),
                   cbind(pri.cons_e_pc,promoter=rep("both"),enhancer=rep("+"),gentype=rep("pc"),gr=rep("pri")), 
                   cbind(pri.cons_ne_pc,promoter=rep("both"),enhancer=rep("-"),gentype=rep("pc"),gr=rep("pri")),
                   cbind(pri.cons_e.np_pc,promoter=rep("-"),enhancer=rep("+"),gentype=rep("pc"),gr=rep("pri")),
                   cbind(pri.cons_ne.p_pc,promoter=rep("+"),enhancer=rep("-"),gentype=rep("pc"),gr=rep("pri")),
                   cbind(pri.cons_e.p_pc,promoter=rep("+"),enhancer=rep("+"),gentype=rep("pc"),gr=rep("pri")),
                   cbind(pri.cons_ne.np_pc,promoter=rep("-"),enhancer=rep("-"),gentype=rep("pc"),gr=rep("pri")))
write.table(x = whole_cons,file = "seq_conserv/enhancer_promoter/whole_cons.txt",quote = F,sep = "\t",row.names = F,col.names = T)

