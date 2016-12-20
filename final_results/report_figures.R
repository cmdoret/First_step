# In this script was used to build figures for the final report
# Cyril Matthey-Doret
# Thu Dec  1 19:14:34 2016 ------------------------------

library(knitr);library(ggplot2)

######################
# Figure: Expression #
######################

#path <- "/Users/cmatthe5/Documents/First_step/data/"
path <- "/home/cyril/Documents/Master/sem_1/First_step/data/"
#path <- "/home/cyril/Documents/First_step/data/"
setwd(path)
cell_lines <- c("GM12878")  # Cell lines to plot
test_lines <- c("GM12878","HUVEC","K562","NHEK")  # Cell lines to include in the table with corr statistics
#test_lines <- c("GM12878")  # Cell lines to include in the table with corr statistics
whole_expr <- read.table("expression/enhancer_promoter/whole_exp.txt",header=T)
whole_expr <- whole_expr[whole_expr$promoter!="both",]
colnames(whole_expr) <- c("gene_id","GM12878","HUVEC","K562","NHEK","promoter","enhancer","gentype")

# Loading lincRNAs sets
e.np_linc <- read.table("enhancer_bound/all_combinations/e.np_linc_pr.bed")  # Overlap enhancer marks, but no promoter marks
ne.np_linc <- read.table("enhancer_bound/all_combinations/ne.np_linc_pr.bed")  # Overlap neither promoter marks, nor enhancer marks

# Loading pcgenes sets
PCG <- read.table("pc_genes/LCL.expressed.pcgene.bed")  # All protein coding genes


#adding colnames
colnames(e.np_linc)=colnames(ne.np_linc)=colnames(PCG) <-c("chr", "start", "end", "gene", "strand")

extract_int <- function(cset, gset){
  csubset <- cset[cset$gene %in% gset$gene,]
  csubset <- csubset[csubset$diam!=0,]
  return(csubset)
}


genlist <- list(e.np_linc,ne.np_linc,PCG)
gennames <- c("elincRNA","other lincRNA","PCG")


corr_expr <- function(hic, expr, cell_line="", gtype=""){  # This function computes the correlation between expression and DNA contacts for a set of genes.
  # Both arguments need to be a list of genes containing the gene ID and its median expression level/amount of DNA-DNA contact.
  merged <- merge(x = hic, y = expr, by = "gene", all = F)
  merged$expression <- as.numeric(merged$expression)
  merged$diam <- as.numeric(merged$diam)
  results <- c(round(cor.test(merged$diam, merged$expression, method="spearman")$estimate,4), 
               round(cor.test(merged$diam, merged$expression, method="spearman")$p.value,4))
  return(results)
}

comp_lines <- data.frame()
# Building a concatenated dataframe containing each gene as many time as there are cell line. More convenient for graphics
for(c in test_lines){
  tmp.expr <- whole_expr[,c("gene_id",c)]
  colnames(tmp.expr) <- c("gene","expression")
  tmp.expr$cell.line <- c
  comp_lines <- rbind(comp_lines,tmp.expr)
}
comp_lines$gentype <- NA  # Preallocating rows for gentypes
count <- 1
for(n in genlist){
  comp_lines[which(comp_lines$gene %in% n$gene),"gentype"] <- gennames[count]
  count <- count+1
}
comp_lines <- na.omit(comp_lines)
# removing rows with empty entries, in case not all gentypes are used (here 3 of 6 are used)
rownames(comp_lines) <- NULL
# resetting row numbers for convenience

# Computing p-values for pairwise comparisons
short_wilcox <- function(x,y){return(format(wilcox.test(x, y)$p.value,digits=3))}
pl<-c()
gl<-c()
tl<-c()
comp_lines <- comp_lines[comp_lines$expression!=0,]
for(t in levels(as.factor(comp_lines$cell.line))){
  p1 <- short_wilcox(comp_lines[comp_lines$cell.line==t & comp_lines$gentype=="elincRNA","expression"],
                     comp_lines[comp_lines$cell.line==t &  comp_lines$gentype=="other lincRNA","expression"])
  p2 <- short_wilcox(comp_lines[comp_lines$cell.line==t & comp_lines$gentype=="elincRNA","expression"],
                     comp_lines[comp_lines$cell.line==t &  comp_lines$gentype=="PCG","expression"])
  pl<-append(pl,c(p1,p2))
  tl<-append(tl,rep(t,2))
  gl<-append(gl,c("elincRNA ~ other lincRNA","elincRNA ~ PCG"))
}
wilcox_p <- data.frame(pval=pl,cell.line=tl,comp=gl)  # Dataframe containing p-values for all pairwise comparisons
starcode <- c()
wilcox_p$pval <- as.character(wilcox_p$pval)
wilcox_p$pval <- as.numeric(wilcox_p$pval)
for(i in wilcox_p$pval){
  print(i)
  if(i<0.1){p <- "."}
  else{p <- "-"}
  if(i<0.05){p <- "*"}
  if(i<0.01){p <- "**"}
  if(i<0.001){p <- "***"}
  print(p)
  starcode <-append(starcode,p)}
wilcox_p$starcode <- starcode
# Building dataframes for displaying bars indicating which boxes correspond to the p-values
bar_pp <- c() # Large bar for elinc_PCG comparison
for(i in seq(0,3*length(test_lines)-1,3)){bar_pp <- append(bar_pp,rep(i,5))}
df1 <- data.frame(a = rep(c(1.1, 1.1,2,2.9,2.9),length(test_lines))+bar_pp, b = c(6.3, 6.5, 6.5, 6.5, 6.3))  
# Path is composed of two small vertical lines and a large horizontal one
bar_p <- c() # short bar for elinc_other linc comparison
for(i in seq(0,3*length(test_lines)-1,3)){bar_p <- append(bar_p,rep(i,4))}
df2 <- data.frame(a = rep(c(1.1, 1.1,1.9, 1.9),length(test_lines))+bar_p, b = c(5.3, 5.5, 5.5, 5.3))


m <- c()
tmp_cell <- c()
tmp_gt <- c()
for(c in test_lines){
  for(g in levels(as.factor(comp_lines$gentype))){
    tmp_cell <- append(tmp_cell, c)
    tmp_gt <- append(tmp_gt, g)
    m <- append(m,median(comp_lines$expression[comp_lines$cell.line==c & comp_lines$gentype==as.character(g)]))
  }
}
med_df <- data.frame(cell.lines = tmp_cell,
                     gentype = tmp_gt,
                     med = m)

increm <- c()
for(i in seq(1,3*length(test_lines),3)){increm <- append(increm,rep(i,2))}  # building vector for positioning p-values between boxes
exp <- ggplot(comp_lines)+geom_boxplot(aes(x=paste(cell.line,gentype,sep="_"),y=log10(expression),fill=gentype), notch=T)+
  scale_x_discrete(labels=rep(gennames,length(test_lines))) +
  scale_fill_manual(values=c("#F8766D","#00BFC4","#C77CFF"))+
  guides(fill=FALSE)+
  theme_bw()+xlab("")+ylab("Log10 expression (RPKM)")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  annotate(geom = 'text',x=rep(c(0.5,1),length(test_lines))+increm, y=rep(c(5.7,6.7),length(test_lines)),label=wilcox_p$starcode,size=5)+
  geom_line(data = df1[1:5,], aes(x = a, y = b)) +
  geom_line(data = df2[1:4,], aes(x = a, y = b)) +
  geom_line(data = df1[6:10,], aes(x = a, y = b)) +
  geom_line(data = df2[5:8,], aes(x = a, y = b)) +
  geom_line(data = df1[11:15,], aes(x = a, y = b)) +
  geom_line(data = df2[9:12,], aes(x = a, y = b)) +
  geom_line(data = df1[16:20,], aes(x = a, y = b)) +
  geom_line(data = df2[13:16,], aes(x = a, y = b))+
  geom_label(data=med_df,x=1:3*length(test_lines),y=log10(as.numeric(med_df$med)),label=round(log10(as.numeric(med_df$med)),3))


########################
# Figure: Conservation #
########################

library(ggplot2);library(gridExtra);library(plyr)

whole_cons <- read.table("../data/seq_conserv/enhancer_promoter/whole_cons.txt", header=T)
whole_cons <- whole_cons[whole_cons$promoter!="both",]
whole_cons2 <- whole_cons[whole_cons$gentype=="pc",]
whole_cons2$gentype <- "PCG"

whole_cons <- whole_cons[whole_cons$promoter=="-" & whole_cons$gentype=="linc",]
ez_class <- c()
rownames(whole_cons) <- NULL
for(r in rownames(whole_cons)){
  r <- as.numeric(r)
  ez_class <- append(ez_class,ifelse(whole_cons$enhancer[r]=="+","elincRNA","other lincRNA"))
}
whole_cons$gentype <- ez_class
whole_cons <- rbind(whole_cons,whole_cons2)
options(digits=3,scipen=0)

#dataf frame for annotation of median
med.fac <- ddply(whole_cons, .(gentype, gr), function(.d)
  data.frame(x=median(na.rm = T,round(.d$avg_score,3))))
#data frame for annotation of p-value
short_wilcox <- function(x,y){return(format(wilcox.test(x, y)$p.value,digits=3))}
pl<-c()
gl<-c()
tl<-c()
for(t in levels(whole_cons$gr)){
  p1 <- short_wilcox(whole_cons[whole_cons$gr==t & whole_cons$gentype=="elincRNA","avg_score"],
                    whole_cons[whole_cons$gr==t &  whole_cons$gentype=="other lincRNA","avg_score"])
  p1 <- as.numeric(p1)
  
  p2 <- short_wilcox(whole_cons[whole_cons$gr==t & whole_cons$gentype=="elincRNA","avg_score"],
                     whole_cons[whole_cons$gr==t &  whole_cons$gentype=="PCG","avg_score"])
  p2 <- as.numeric(p2)
  if(nchar(p1)>4){p1 <-format(p1,scientific=T)}
  if(nchar(p2)>4){p2 <-format(p2,scientific=T)}
  pl<-append(pl,c(p1,p2))
  tl<-append(tl,rep(t,2))
  gl<-append(gl,c("elincRNA ~ other lincRNA","elincRNA ~ PCG"))
}
wilcox_p <- data.frame(pval=pl,gr=tl,gentype=gl)
starcode <- c()
wilcox_p$pval <- as.character(wilcox_p$pval)
wilcox_p$pval <- as.numeric(wilcox_p$pval)
for(i in wilcox_p$pval){
  if(i<0.1){p <- "."}
  else{p <- "-"}
  if(i<0.05){p <- "*"}
  if(i<0.01){p <- "**"}
  if(i<0.001){p <- "***"}
  starcode <-append(starcode,p)}
wilcox_p$starcode <- starcode
# Building dataframes for displaying bars indicating which boxes correspond to the p-values
bar_pp <- c() # Large bar for elinc_PCG comparison
for(i in seq(0,3*length(test_lines)-1,3)){bar_pp <- append(bar_pp,rep(i,5))}
df1 <- data.frame(a = rep(c(1.1, 1.1,2,2.9,2.9),length(test_lines))+bar_pp, b = c(0.97, 1, 1, 1, 0.97))  
# Path is composed of two small vertical lines and a large horizontal one
bar_p <- c() # short bar for elinc_other linc comparison
for(i in seq(0,3*length(test_lines)-1,3)){bar_p <- append(bar_p,rep(i,4))}
df2 <- data.frame(a = rep(c(1.1, 1.1,1.9, 1.9),length(test_lines))+bar_p, b = c(0.9, 0.93, 0.93, 0.9))



arcons <- read.table("seq_conserv/enhancer_bound/whole_cons.txt",header=T)
arcons <- arcons[arcons$gentype=="AR",]
group_names <- c(`mam` = "Mammals",`pri` = "Primates")
gentype_names <- c(`lincRNA` = "lincRNAs",`pc` = "protein-coding")
AR_medcons <- data.frame(val=c(median(arcons$avg_score[arcons$gentype=="AR" & arcons$gr=="mam"],na.rm = T),
                               median(arcons$avg_score[arcons$gentype=="AR" & arcons$gr=="pri"],na.rm = T)),
                         gr=c("mam","pri"))

increm <- c()
for(i in seq(1,3*length(test_lines),3)){increm <- append(increm,rep(i,2))}  # building vector for positioning p-values between boxes
cons <- ggplot(data=whole_cons)+
  facet_grid(~gr,labeller = labeller(gr = as_labeller(group_names)))+
  geom_boxplot(aes(x=gentype,y=avg_score,fill=gentype),notch=T)+
  geom_hline(data= AR_medcons,aes(yintercept = AR_medcons$val),size=1,col="#66bb66",show.legend= F)+
  scale_fill_manual(values = c("#F8766D","#00BFC4","#C77CFF"))+
  #scale_linetype_manual("Title", values = 2) +
  #guides(fill = guide_legend(title = element_blank(),override.aes=list(linetype=0)))+
  #geom_text(data=wilcox_p, aes(x=c(1.5, 2), y=1, label=paste0("p=",pval)), 
  #          colour="black", inherit.aes=T, parse=T)+
  geom_label(data=med.fac, aes(x=gentype, y=x+0.01, label=x), 
             colour="black", inherit.aes=FALSE, parse=FALSE,size=3)+
  theme_bw()+ ylab("averaged phastCons")+xlab("")+guides(fill=F,line=F)+coord_cartesian(ylim=c(0,1.1))+theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_text(data=wilcox_p,aes(x=rep(c(0.5,1),length(test_lines)*2)+increm, y=rep(c(0.94,1.01),length(test_lines)*2),label=starcode,size=5))+
  geom_line(data = df1[1:5,], aes(x = a, y = b)) +
  geom_line(data = df2[1:4,], aes(x = a, y = b)) +
  geom_line(data = df1[6:10,], aes(x = a, y = b)) +
  geom_line(data = df2[5:8,], aes(x = a, y = b)) +
  geom_line(data = df1[11:15,], aes(x = a, y = b)) +
  geom_line(data = df2[9:12,], aes(x = a, y = b)) +
  geom_line(data = df1[16:20,], aes(x = a, y = b)) +
  geom_line(data = df2[13:16,], aes(x = a, y = b))
grid.arrange(layout_matrix=matrix(c(1,1,1,2,2,2,2,2,2),nrow=3,byrow=F),grobs=list(exp,cons),nrow=1)
##############################
# Figure: Tissue specificity #
##############################

whole_tau <- read.table("../data/tissue_specificity/enhancer_promoter/whole_tau.txt", header=T)
whole_tau <- whole_tau[whole_tau$promoter!="both",]
whole_tau2 <- whole_tau[whole_tau$gentype=="pc",]
whole_tau2$gentype <- "PCG"

whole_tau <- whole_tau[whole_tau$promoter=="-" & whole_tau$gentype=="linc",]
ez_class <- c()
rownames(whole_tau) <- NULL
for(r in rownames(whole_tau)){
  r <- as.numeric(r)
  ez_class <- append(ez_class,ifelse(whole_tau$enhancer[r]=="+","elincRNA","other lincRNA"))
}
whole_tau$gentype <- ez_class
whole_tau <- rbind(whole_tau,whole_tau2)

options(digits=3,scipen=0)

#dataf frame for annotation of median
library(plyr)
med.fac <- ddply(whole_tau, .(gentype), function(.d)
  data.frame(x=median(na.rm = T,round(.d$tau,3))))
#data frame for annotation of p-value
short_wilcox <- function(x,y){return(format(wilcox.test(x, y)$p.value,digits=3))}
pl<-c()
gl<-c()

p1 <- short_wilcox(whole_tau[whole_tau$gentype=="elincRNA","tau"],
                   whole_tau[whole_tau$gentype=="other lincRNA","tau"])
p1 <- as.numeric(p1)

p2 <- short_wilcox(whole_tau[whole_tau$gentype=="elincRNA","tau"],
                   whole_tau[whole_tau$gentype=="PCG","tau"])
p2 <- as.numeric(p2)
if(nchar(p1)>4){p1 <-format(p1,scientific=T)}
if(nchar(p2)>4){p2 <-format(p2,scientific=T)}
pl<-append(pl,c(p1,p2))
gl<-append(gl,c("elincRNA ~ other lincRNA","elincRNA ~ PCG"))

wilcox_p <- data.frame(pval=pl,gentype=gl)
starcode <- c()
wilcox_p$pval <- as.character(wilcox_p$pval)
wilcox_p$pval <- as.numeric(wilcox_p$pval)
for(i in wilcox_p$pval){
  if(i<0.1){p <- "."}
  else{p <- "-"}
  if(i<0.05){p <- "*"}
  if(i<0.01){p <- "**"}
  if(i<0.001){p <- "***"}
  starcode <-append(starcode,p)}
wilcox_p$starcode <- starcode
# Building dataframes for displaying bars indicating which boxes correspond to the p-values
bar_pp <- c() # Large bar for elinc_PCG comparison
for(i in seq(0,3*length(test_lines)-1,3)){bar_pp <- append(bar_pp,rep(i,5))}
df1 <- data.frame(a = rep(c(1.1, 1.1,2,2.9,2.9),length(test_lines))+bar_pp, b = c(1.07, 1.1, 1.1, 1.1, 1.07))  
# Path is composed of two small vertical lines and a large horizontal one
bar_p <- c() # short bar for elinc_other linc comparison
for(i in seq(0,3*length(test_lines)-1,3)){bar_p <- append(bar_p,rep(i,4))}
df2 <- data.frame(a = rep(c(1.1, 1.1,1.9, 1.9),length(test_lines))+bar_p, b = c(1, 1.03, 1.03, 1))





increm <- c()
for(i in seq(1,3*length(test_lines),3)){increm <- append(increm,rep(i,2))}  # building vector for positioning p-values between boxes
spec <-ggplot(data=whole_tau)+
  geom_boxplot(aes(x=gentype,y=tau,fill=gentype),notch=T)+
  scale_fill_manual(values = c("#F8766D","#00BFC4","#C77CFF"))+
  #scale_linetype_manual("Title", values = 2) +
  #guides(fill = guide_legend(title = element_blank(),override.aes=list(linetype=0)))+
  #geom_text(data=wilcox_p, aes(x=c(1.5, 2), y=1, label=paste0("p=",pval)), 
  #          colour="black", inherit.aes=T, parse=T)+
  geom_label(data=med.fac, aes(x=gentype, y=x+0.01, label=x), 
             colour="black", inherit.aes=FALSE, parse=FALSE,size=3)+
  theme_bw()+ ylab("Tau")+xlab("")+guides(fill=F,line=F)+coord_cartesian(ylim=c(0,1.1))+theme(legend.position="none",axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_text(data=wilcox_p,aes(x=rep(c(0.5,1),length(test_lines))+increm, y=rep(c(1.05,1.11),length(test_lines)),label=starcode,size=5))+
  geom_line(data = df1[1:5,], aes(x = a, y = b)) +
  geom_line(data = df2[1:4,], aes(x = a, y = b)) +
  geom_line(data = df1[6:10,], aes(x = a, y = b)) +
  geom_line(data = df2[5:8,], aes(x = a, y = b)) +
  geom_line(data = df1[11:15,], aes(x = a, y = b)) +
  geom_line(data = df2[9:12,], aes(x = a, y = b)) +
  geom_line(data = df1[16:20,], aes(x = a, y = b)) +
  geom_line(data = df2[13:16,], aes(x = a, y = b))

grid.arrange(layout_matrix=matrix(c(1,2,3,2),nrow=2,byrow=F),grobs=list(exp,cons,spec),nrow=1)
#####################################################
# Figure: Enrichment at HiC-boundaries and anchors #
#####################################################
library(ggplot2)
whole_hicbound <- read.table("GAT/out/whole_seg_10kgat_hic_boundaries.txt",header=T)
whole_hicbound <- whole_hicbound[whole_hicbound$segment %in% c("e.np", "ne.np"),]
whole_anchors <- read.table("GAT/out/whole_seg_10kgat_anchors.txt",header=T)
whole_anchors <- whole_anchors[whole_anchors$segment %in% c("e.np", "ne.np"),]
nicenames <- c("elincRNA \nprom","elincRNA \nprom+body","other lincRNA \nprom", "other lincRNA \nprom+body")
barplot(whole_hicbound$fold,names.arg = nicenames,main=" Enrichment at Hi-C boundaries")
abline(h = whole_hicbound$fold[1])
abline(h = whole_hicbound$fold[2])
text(x = c(0.5,2,3,4.5),y=whole_hicbound$fold-0.1,labels = whole_hicbound$fold)

library(gridExtra)
g1 <- ggplot(data=whole_hicbound[whole_hicbound$element=="pr",],aes(x=nicenames[c(1,3)],y=fold,fill=segment))+
  geom_bar(stat = 'identity')+
  #scale_fill_continuous(high="#DD5555",low = "#55DD55")+
  theme_bw()+guides(fill=F)+ylab("Fold enrichment")+
  #geom_hline(data=whole_hicbound[1:2,],aes(yintercept=fold[c(1,2)]))+
  xlab("")+ggtitle("Enrichment at TAD boundaries")+
  geom_text(aes(x=c(1,2),y=fold-0.1,label=fold))+
  geom_text(aes(x=c(1,2),y=fold-0.2,label=paste0("q=",qval)))+
  scale_x_discrete(labels=nicenames[c(1,3)])+coord_cartesian(ylim=c(0,3))

g2 <- ggplot(data=whole_anchors[whole_anchors$element=="pr",],aes(x=nicenames[c(1,3)],y=fold,fill=segment))+
  geom_bar(stat = 'identity')+
  #scale_fill_continuous(high="#DD5555",low = "#55DD55")+
  theme_bw()+ylab("Fold enrichment")+guides(fill=F)+
  #geom_hline(data=whole_anchors[1:2,],aes(yintercept=fold[c(1,2)]))+
  xlab("")+ggtitle("Enrichment at loops anchors")+
  geom_text(aes(x=c(1,2),y=fold-0.1,label=fold))+
  geom_text(aes(x=c(1,2),y=fold-0.2,label=paste0("q=",qval)))+
  scale_x_discrete(labels=nicenames[c(1,3)])+coord_cartesian(ylim=c(0,3))
grid.arrange(g1,g2,nrow=1)

whole_bins <- read.table("GAT/out/whole_seg_10kgat_results.txt",header=T)




###################################
# Figure: Enrichment per TAD bins #
###################################
library(gridExtra)
setwd("/home/cyril/Documents/First_step/data/")
#setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
#setwd("/Users/cmatthe5/Documents/First_step/data/")
whole_gat_seg <- read.table("GAT/out/whole_seg_10kgat_fullover_results.txt",header=T)
whole_gat_seg <- whole_gat_seg[whole_gat_seg$element=="pr" &
                                 !(whole_gat_seg$annotation %in% c("ne","e"))&
                                 !(whole_gat_seg$segment %in% c("ne","e")),]
sets <- c("e.np","ne.np")
nice_workspace <- c(wholegenome= "Whole genome",intergenic="Intergenic space",allpc="Protein-coding space", lincRNA="Expressed lincRNAs space")
rownames(whole_gat_seg) <- NULL
gene2bins_seg<-droplevels(whole_gat_seg[whole_gat_seg$annotation=="short10bins" &
                                          (whole_gat_seg$segment %in% sets | whole_gat_seg$workspace=="allpc"),]) 
gene2bins_seg$annottrack <- as.character(gene2bins_seg$annottrack)
gene2bins_seg$annottrack <- factor(gene2bins_seg$annottrack,levels = c(1:20))

binplot<-ggplot(data=gene2bins_seg[gene2bins_seg$workspace=="intergenic",],
       aes(x=annottrack,y=fold,fill=log10(qval)))+
  facet_grid(~segment, labeller = as_labeller(c(`e.np`="elincRNA",`ne.np`="other lincRNA")))+
  geom_bar(stat = 'identity')+
  scale_fill_continuous(high="#DD5555",low = "#55DD55",name="Log10 q-value",
                        limits = c(-3,0), breaks=c(0,-1, -2, -3),guide = guide_colorbar(nbin=5))+
  theme_bw()+ylab("Fold enrichment")+
  #geom_hline(data=whole_anchors[1:2,],aes(yintercept=fold[c(1,2)]))+
  xlab("")+ggtitle("Enrichment of lincRNAs across TAD bins")
  #geom_text(aes(x=c(1,2),y=fold-0.3,label=fold))+
  #geom_text(aes(x=c(1,2),y=fold-0.8,label=paste0("q=",qval)))+



################################################
# Figure: Enrichment in architectural proteins #
################################################

chip_gene_seg<-droplevels(whole_gat_seg[(whole_gat_seg$annotation %in% sets | whole_gat_seg$workspace=="allpc") &
                                           whole_gat_seg$segment %in% c("GM12878RAD21","GM12878SMC3","GM12878CTCF"),])

ggplot(data=chip_gene_seg[chip_gene_seg$workspace=="intergenic",],aes(x=annotation,y=fold,fill=annotation))+
  geom_bar(stat = 'identity')+
  facet_grid(~segment,labeller = as_labeller(c(`GM12878RAD21`="RAD21",`GM12878SMC3`="SMC3",`GM12878CTCF`="CTCF")))+
  #scale_fill_continuous(high="#DD5555",low = "#55DD55")+
  theme_bw()+ylab("Fold enrichment")+guides(fill=F)+
  #geom_hline(data=whole_anchors[1:2,],aes(yintercept=fold[c(1,2)]))+
  xlab("")+ggtitle("Enrichment of architectural proteins\n binding sites")+
  geom_text(aes(y=fold-0.3,label=fold))+
  geom_text(aes(y=fold-0.8,label=paste0("q=",qval)))+
  scale_x_discrete(labels=nicenames[c(1,3)])+coord_cartesian(ylim=c(0,8))


library(ggplot2);library(gridExtra)


###########################################################
# Figure: Enrichment in mutually exclusive binding sites #
###########################################################
library(gridExtra)
whole_exclu <- read.table("GAT/out/whole_seg_10kgat_exclu_insul_bs.txt",header=T)
whole_exclu <- whole_exclu[whole_exclu$track!="merged",]
whole_exclu <- whole_exclu[whole_exclu$annotation %in% c("e.np", "ne.np"),]

e1 <-ggplot(data=whole_exclu[whole_exclu$segment=="exclu.CTCF" & whole_exclu$element=="pr",],
       aes(x=nicenames[c(1,3)],y=fold,fill=annotation))+
  geom_bar(stat = 'identity')+
  #scale_fill_continuous(high="#DD5555",low = "#55DD55")+
  theme_bw()+ylab("Fold enrichment")+guides(fill=F)+
  #geom_hline(data=whole_anchors[1:2,],aes(yintercept=fold[c(1,2)]))+
  xlab("")+ggtitle("Enrichment of CTCF-only\n binding sites")+
  geom_text(aes(x=c(1,2),y=fold-0.3,label=fold))+
  geom_text(aes(x=c(1,2),y=fold-0.8,label=paste0("q=",qval)))+
  scale_x_discrete(labels=nicenames[c(1,3)])+coord_cartesian(ylim=c(0,14))

e2 <-ggplot(data=whole_exclu[whole_exclu$segment=="exclu.cohesin" & whole_exclu$element=="pr",],
       aes(x=nicenames[c(1,3)],y=fold,fill=annotation))+
  geom_bar(stat = 'identity')+
  #scale_fill_continuous(high="#DD5555",low = "#55DD55")+
  theme_bw()+ylab("Fold enrichment")+guides(fill=F)+
  #geom_hline(data=whole_anchors[1:2,],aes(yintercept=fold[c(1,2)]))+
  xlab("")+ggtitle("Enrichment of cohesin-only\n binding sites")+
  geom_text(aes(x=c(1,2),y=fold-0.3,label=fold))+
  geom_text(aes(x=c(1,2),y=fold-0.8,label=paste0("q=",qval)))+
  scale_x_discrete(labels=nicenames[c(1,3)])+coord_cartesian(ylim=c(0,14))
grid.arrange(e1,e2,nrow=1)

excluplot <-ggplot(data=whole_exclu[whole_exclu$element=="pr",],aes(x=nicenames[c(1,3,1,3)],y=fold,fill=annotation))+
  geom_bar(stat = 'identity')+
  facet_grid(~segment,labeller = as_labeller(c(`exclu.cohesin`="Cohesin-only",`exclu.CTCF`="CTCF-only")))+
  #scale_fill_continuous(high="#DD5555",low = "#55DD55")+
  theme_bw()+ylab("Fold enrichment")+guides(fill=F)+
  #geom_hline(data=whole_anchors[1:2,],aes(yintercept=fold[c(1,2)]))+
  xlab("")+ggtitle("Enrichment of architectural proteins exclusive binding sites")+
  geom_text(aes(y=fold-0.3,label=fold))+
  geom_text(aes(y=fold-0.8,label=paste0("q=",qval)))+
  scale_x_discrete(labels=nicenames[c(1,3)])+coord_cartesian(ylim=c(0,14))


#############################################################
# Figure: Overlap between SMC3, RAD21 and CTCF binding sites#
#############################################################

setwd("/Users/cmatthe5/Documents/First_step/data/")

# Loading bindingsites and overlaps
RAD21 <- read.table("chip_seq/GM12878_RAD21_peaks.bed")
CTCF <- read.table("chip_seq/GM12878_CTCF_peaks.bed")
SMC3 <- read.table("chip_seq/GM12878_SMC3_peaks.bed")
RAD21_CTCF <- read.table("chip_seq/inter_insulators_GM12878/inter_CTCF_RAD21.bed")
RAD21_SMC3 <- read.table("chip_seq/inter_insulators_GM12878/inter_RAD21_SMC3.bed")
SMC3_CTCF <- read.table("chip_seq/inter_insulators_GM12878/inter_CTCF_SMC3.bed")
all_inter <- read.table("chip_seq/inter_insulators_GM12878/inter_CTCF_RAD21_SMC3.bed")

# Building Venn diagram:
library(VennDiagram)

grid.newpage()
vennplot<-draw.triple.venn(euler.d = T,scaled = T, area1 = length(RAD21$V1), area2 = length(CTCF$V1), area3 = length(SMC3$V1), n12 = length(RAD21_CTCF$V1), 
                 n23 = length(SMC3_CTCF$V1), n13 = length(RAD21_SMC3$V1), 
                 n123 = length(all_inter$V1), category = c("RAD21", "CTCF", "SMC3"), lty = "blank", 
                 fill = c("skyblue", "yellow", "mediumorchid"))

#######################
# Figure: HiC contact #
#######################
library(ggplot2)
setwd(path)

ref_set <- read.table("../data/linc_RNA/LCL.expressed.lincRNA.bed")
com.gen <- rep(NA,length(ref_set[,1]))
comp_lines <- data.frame()
# Building a concatenated dataframe containing each gene as many time as there are cell line. More convenient for graphics
for(c in test_lines){
  tmp.contact <- read.table(paste0("TAD_contact/TAD_TAD_contact/all.lincRNA.",c,".HiC.contact.txt"), header=T)
  tmp.contact$cell.line <- c
  comp_lines <- rbind(comp_lines,tmp.contact)
}
comp_lines$gentype <- NA  # Preallocating rows for gentypes
ne.np_linc <- read.table("enhancer_bound/all_combinations/ne.np_linc_pr.bed")
e.np_linc <- read.table("enhancer_bound/all_combinations/e.np_linc_pr.bed")
colnames(e.np_linc) = colnames(ne.np_linc)<-c("chr", "start", "end", "gene", "strand")
linclist <- list(e.np_linc,ne.np_linc)
lincnames <- c("elincRNA","other lincRNA")
count <- 1
for(n in linclist){
  comp_lines[which(comp_lines$gene %in% n$gene),"gentype"] <- lincnames[count]
  count <- count+1
}
comp_lines <- na.omit(comp_lines)
# removing rows with empty entries, in case not all gentypes are used (here 3 of 6 are used)
rownames(comp_lines) <- NULL
# resetting row numbers for convenience

# Computing p-values for pairwise comparisons
short_wilcox <- function(x,y){return(format(wilcox.test(x, y)$p.value,digits=3))}
pl<-c()
gl<-c()
tl<-c()
comp_lines <- comp_lines[comp_lines$diam!=0,]
for(t in levels(as.factor(comp_lines$cell.line))){
  p1 <- short_wilcox(comp_lines[comp_lines$cell.line==t & comp_lines$gentype=="elincRNA","diam"],
                     comp_lines[comp_lines$cell.line==t &  comp_lines$gentype=="other lincRNA","diam"])
#  p2 <- short_wilcox(comp_lines[comp_lines$cell.line==t & comp_lines$gentype=="enplincRNA","diam"],
#                     comp_lines[comp_lines$cell.line==t &  comp_lines$gentype=="neplincRNA","diam"])
  pl<-append(pl,p1)
  tl<-append(tl,t)
  gl<-append(gl,"elincRNA ~ other lincRNA")
}
wilcox_p <- data.frame(pval=pl,cell.line=tl,comp=gl)  # Dataframe containing p-values for all pairwise comparisons
starcode <- c()
wilcox_p$pval <- as.character(wilcox_p$pval)
wilcox_p$pval <- as.numeric(wilcox_p$pval)
for(i in wilcox_p$pval){
  if(i<0.1){p <- "."}
  else{p <- "-"}
  if(i<0.05){p <- "*"}
  if(i<0.01){p <- "**"}
  if(i<0.001){p <- "***"}
  starcode <-append(starcode,p)}
wilcox_p$starcode <- starcode
# Building dataframes for displaying bars indicating which boxes correspond to the p-values
#bar_pp <- c()
#for(i in seq(0,11,3)){bar_pp <- append(bar_pp,rep(i,5))}
#df1 <- data.frame(a = rep(c(1.1, 1.1,2,2.9,2.9),4)+bar_pp, b = c(73, 75, 75, 75, 73))
bar_p <- c()
for(i in seq(0,7,2)){bar_p <- append(bar_p,rep(i,4))}
df2 <- data.frame(a = rep(c(1.1, 1.1,1.9, 1.9),4)+bar_p, b = c(63, 65, 65, 63))
#df3 <- data.frame(a = rep(c(2.1, 2.1, 2.9, 2.9),4)+bar_p, b = c(53, 55, 55, 53))

increm <- c()
for(i in seq(1,8,2)){increm <- append(increm,rep(i,1))}  # building vector for positioning p-values between boxes
ggplot(comp_lines)+geom_boxplot(aes(x=paste(cell.line,gentype,sep="_"),y=diam,fill=cell.line))+
  scale_x_discrete(labels=rep(lincnames[c(1,2)],4)) +
  theme_bw()+xlab("Gene class")+ylab("TAD vs TAD contacts")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  annotate(geom = 'text',x=rep(c(0.5),4)+increm, y=rep(c(67),4),label=wilcox_p$starcode,size=5)+
  #geom_line(data = df1[1:5,], aes(x = a, y = b)) +
  geom_line(data = df2[1:4,], aes(x = a, y = b)) +
  #geom_line(data = df3[1:4,], aes(x = a, y = b)) +
  #geom_line(data = df1[6:10,], aes(x = a, y = b)) +
  geom_line(data = df2[5:8,], aes(x = a, y = b)) +
  #geom_line(data = df3[5:8,], aes(x = a, y = b)) +
  #geom_line(data = df1[11:15,], aes(x = a, y = b)) +
  geom_line(data = df2[9:12,], aes(x = a, y = b)) +
  #geom_line(data = df3[9:12,], aes(x = a, y = b)) +
  #geom_line(data = df1[16:20,], aes(x = a, y = b)) +
  geom_line(data = df2[13:16,], aes(x = a, y = b))
  #geom_line(data = df3[13:16,], aes(x = a, y = b))

####################
# Presentation only: simplified enrichment: RAD21 and SMC3 merged
####################
sets <- c("e.np","ne.np")
union_gat <- read.table("GAT/out/whole_seg_10kgat_cohesin_union_results.txt",header=T, sep="\t")
union_gat<-droplevels(union_gat[(union_gat$annotation %in% sets | union_gat$workspace=="allpc") &
                                          union_gat$segment == "unionSMC3RAD21merged",])
union_gat <- union_gat[union_gat$element=="pr",]
union_gat <- rbind(union_gat, chip_gene_seg[chip_gene_seg$segment=="GM12878CTCF" & chip_gene_seg$workspace=="intergenic",])
union <- ggplot(data=union_gat,aes(x=annotation,y=fold,fill=annotation))+
  geom_bar(stat = 'identity')+
  facet_grid(~segment,labeller = as_labeller(c(`unionSMC3RAD21merged`="Cohesin",`GM12878CTCF`="CTCF")))+
  #scale_fill_continuous(high="#DD5555",low = "#55DD55")+
  theme_bw()+ylab("Fold enrichment")+guides(fill=F)+
  #geom_hline(data=whole_anchors[1:2,],aes(yintercept=fold[c(1,2)]))+
  xlab("")+ggtitle("Enrichment of architectural proteins binding sites")+
  geom_text(aes(y=fold-0.3,label=fold))+
  geom_text(aes(y=fold-0.8,label=paste0("q=",qval)))+
  scale_x_discrete(labels=nicenames[c(1,3)])+coord_cartesian(ylim=c(0,14))

grid.arrange(layout_matrix=matrix(c(1,2,3,4),ncol=2,byrow=F),grobs=list(anchorplot,union,gTree(children=vennplot),excluplot))
grid.arrange(layout_matrix=matrix(c(1,2),ncol=2,byrow=T),grobs=list(boundplot,binplot))
library(ggplot2);library(gridExtra)
