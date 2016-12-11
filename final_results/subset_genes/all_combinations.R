# This script is used after bedtools intersect and enhanc_assoc calls to classify genes depending on if they overlap promoters AND 
# enhancers (all combinations possible) as defined in the report.
# Cyril Matthey-Doret
# Fri Nov 18 14:14:27 2016 ------------------------------


setwd("/home/cyril/Documents/Master/sem_1/First_step/data/enhancer_bound/all_combinations/")

nelinc <- read.table("ne_linc_pr.bed")  # lincRNAs not overlapping enhancers
nepc <- read.table("ne_pc_pr.bed")  # PCG not overlapping enhancers
pnelinc <-read.table("ne.p_linc_pr.bed")  # lincRNAs overlapping promoters but no enhancers
pnepc <- read.table("ne.p_pc_pr.bed")  # PCG overlapping promoter but no enhancers
colnames(nepc)=colnames(nelinc)=colnames(pnepc)=
  colnames(pnelinc) <- c("chr", "start", "end","gene","strand")

npne_linc <- nelinc[!(nelinc$gene %in% pnelinc$gene),] # lincRNAs overlapping neither enhancers nor promoters
npne_pc <- nepc[!(nepc$gene %in% pnepc$gene),]  # PCG overlapping neither enhancers nor promoters

write.table(npne_linc,"ne.np_linc_pr.bed",col.names=F,row.names = F,quote=F,sep="\t") # Other lincRNAs
write.table(npne_pc,"ne.np_pc_pr.bed",col.names=F,row.names = F,quote=F,sep="\t")


elinc <- read.table("e_linc_pr.bed")  # lincRNAs overlapping enhancers
epc <- read.table("e_pc_pr.bed")  # PCG overlapping enhancers
pelinc <-read.table("e.p_linc_pr.bed")  # lincRNAs overlapping both promoters and enhancers
pepc <- read.table("e.p_pc_pr.bed")  # PCG overlapping both promoters and enhancers
colnames(epc)=colnames(elinc)=colnames(pepc)=
  colnames(pelinc) <- c("chr", "start", "end","gene","strand")

npe_linc <- elinc[!(elinc$gene %in% pelinc$gene),]  # lincRNAs overlapping enhancers but no promoters
npe_pc <- epc[!(epc$gene %in% pepc$gene),]  # PVG overlapping enhancers but no promoters

write.table(npe_linc,"e.np_linc_pr.bed",col.names=F,row.names = F,quote=F,sep="\t")  # ElincRNAs
write.table(npe_pc,"e.np_pc_pr.bed",col.names=F,row.names = F,quote=F,sep="\t")

