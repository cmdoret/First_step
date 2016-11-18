# Here I categorize lincRNAs and protein-coding genes depending on if they overlap promoters and enhancers (all combinations possible)
# Cyril
# Fri Nov 18 14:14:27 2016 ------------------------------


setwd("/home/cyril/Documents/Master/sem_1/First_step/data/enhancer_bound/all_combinations/")

nelinc <- read.table("ne_linc_pr.bed")
nepc <- read.table("ne_pc_pr.bed")
pnelinc <-read.table("ne.p_linc_pr.bed")
pnepc <- read.table("ne.p_pc_pr.bed")
colnames(nepc)=colnames(nelinc)=colnames(pnepc)=
  colnames(pnelinc) <- c("chr", "start", "end","gene","strand")

npne_linc <- nelinc[!(nelinc$gene %in% pnelinc$gene),]
npne_pc <- nepc[!(nepc$gene %in% pnepc$gene),]

write.table(npne_linc,"ne.np_linc_pr.bed",col.names=F,row.names = F,quote=F,sep="\t")
write.table(npne_pc,"ne.np_pc_pr.bed",col.names=F,row.names = F,quote=F,sep="\t")


elinc <- read.table("e_linc_pr.bed")
epc <- read.table("e_pc_pr.bed")
pelinc <-read.table("e.p_linc_pr.bed")
pepc <- read.table("e.p_pc_pr.bed")
colnames(epc)=colnames(elinc)=colnames(pepc)=
  colnames(pelinc) <- c("chr", "start", "end","gene","strand")

npe_linc <- elinc[!(elinc$gene %in% pelinc$gene),]
npe_pc <- epc[!(epc$gene %in% pepc$gene),]

write.table(npe_linc,"e.np_linc_pr.bed",col.names=F,row.names = F,quote=F,sep="\t")
write.table(npe_pc,"e.np_pc_pr.bed",col.names=F,row.names = F,quote=F,sep="\t")

