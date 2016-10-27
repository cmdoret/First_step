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
# pr = promoter region; prb = promoter region + body
gene_2_pr <- function(g,b=FALSE){
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

stupid_loop <-function(g,b=FALSE){
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

