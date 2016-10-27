# This script splits lincRNAs and pcgenes into enhancer associating and non-enhancer associated.
# Cyril Matthey-Doret
# Thu Oct 27 15:40:10 2016 ------------------------------

# Loading data:
setwd("/Users/cmatthe5/Documents/First_step/data/")
linc <- read.table("linc_RNA/LCL.expressed.lincRNA.bed")
pc <- read.table("pc_genes/LCL.expressed.pcgene.bed")
colnames(pc) = colnames(linc) <- c("chr","start","end","gene","strand")

#=============================
# Editing intervals for inclusion of inclusion of promoter region
# Functions defined below allow to extend all genes by 1kb on the side of the promoter
# or to select only the promoter region.
# pr = promoter region; prb = promoter region + body
gene_2_pr <- function(g,b=FALSE){  # This function breaks the bed files DO NOT USE!
  ps <- 0
  if(g[5]=="-"){
    ps <- g[3]
    g[3] <- as.numeric(g[3])+1000
    if(b==F){g[2] <- ps}
  }else{
    ps <- as.numeric(g[2])-1000
    if(b==F){g[3] <- g[2]}
    g[2] <- ps
  }
  return(g)
}

stupid_loop <-function(g,b=FALSE){ # This one is fine
  for(r in 1:length(g[,1])){
    if(g[r,5]=="-"){
      ps <- g[r,3]
      g[r,3] <- as.numeric(g[r,3]+1000)
      if(b==FALSE){g[r,2] <- ps}
    } else{
      ps <- as.numeric(g[r,2])-1000
      if(b==F){g[r,3] <- g[r,2]}
      g[r,2] <- ps
    }
  }
  return(g)
}

#=============================
# Defining promoter region and promoter region + body of lincRNA, pcgenes

pr_linc <- t(apply(X=linc,MARGIN = 1,FUN=gene_2_pr))
prb_linc <- t(apply(X=linc,MARGIN = 1,FUN=gene_2_pr,b=TRUE))
pr_pc <- t(apply(X=pc,MARGIN = 1,FUN=gene_2_pr))
prb_pc <- t(apply(X=pc,MARGIN = 1,FUN=gene_2_pr,b=TRUE))
# This doesn't work, probably this sneaky apply that changes smth in the file...
#============================
# Stupid loop version
pr_linc <- stupid_loop(linc)
prb_linc <- stupid_loop(linc,b=TRUE)
pr_pc <- stupid_loop(pc)
prb_pc <- stupid_loop(pc,b=TRUE)


#============================
# Writing into files
write.table(pr_linc,file="enhancer_bound/pro_linc.bed",quote=F,sep="\t",col.names = F,row.names = F)
write.table(prb_linc,file="enhancer_bound/probo_linc.bed",quote=F,sep="\t",col.names = F,row.names = F)
write.table(pr_pc,file="enhancer_bound/pro_pc.bed",quote=F,sep="\t",col.names = F,row.names = F)
write.table(prb_pc,file="enhancer_bound/probo_pc.bed",quote=F,sep="\t",col.names = F,row.names = F)

#============================
# BEDtools intersect with ENCODE enhancers elements...  

#============================
# Loading overlaps:

over_pr_linc <- read.table("enhancer_bound/overlaps/over_enhancer_pro_linc.bed")
over_prb_linc <- read.table("enhancer_bound/overlaps/over_enhancer_probo_linc.bed")
over_pr_pc <- read.table("enhancer_bound/overlaps/over_enhancer_pro_pc.bed")
over_prb_pc <- read.table("enhancer_bound/overlaps/over_enhancer_probo_pc.bed")
colnames(over_pr_pc)=colnames(over_prb_pc)=colnames(over_pr_linc)=
  colnames(over_prb_linc) <- c("chr","start","end","gene","strand")
#===========================
# Removing duplicates to identify enhancer-associated genes:
# abbreviations: e=enhancer-associated, ne=non-enhancer-associated

elinc_pr <- over_pr_linc[!duplicated(over_pr_linc$gene),]
elinc_prb <- over_prb_linc[!duplicated(over_prb_linc$gene),]
epc_pr <- over_pr_pc[!duplicated(over_pr_pc$gene),]
epc_prb <- over_prb_pc[!duplicated(over_prb_pc$gene),]
write.table(elinc_pr,file = "enhancer_bound/elinc_pr.bed",sep="\t",quote=F,col.names=F,row.names = F)
write.table(elinc_prb,file = "enhancer_bound/elinc_prb.bed",sep="\t",quote=F,col.names=F,row.names = F)
write.table(epc_pr,file = "enhancer_bound/epc_pr.bed",sep="\t",quote=F,col.names=F,row.names = F)
write.table(epc_prb,file = "enhancer_bound/epc_prb.bed",sep="\t",quote=F,col.names=F,row.names = F)

nelinc_pr <- pr_linc[!(pr_linc$gene %in% elinc_pr$gene),]
nelinc_prb <- prb_linc[!(prb_linc$gene %in% elinc_prb$gene),]
nepc_pr <- pr_pc[!(pr_pc$gene %in% epc_pr$gene),]
nepc_prb <- prb_pc[!(prb_pc$gene %in% epc_prb$gene),]

write.table(nelinc_pr,file = "enhancer_bound/nelinc_pr.bed",sep="\t",quote=F,col.names=F,row.names = F)
write.table(nelinc_prb,file = "enhancer_bound/nelinc_prb.bed",sep="\t",quote=F,col.names=F,row.names = F)
write.table(nepc_pr,file = "enhancer_bound/nepc_pr.bed",sep="\t",quote=F,col.names=F,row.names = F)
write.table(nepc_prb,file = "enhancer_bound/nepc_prb.bed",sep="\t",quote=F,col.names=F,row.names = F)