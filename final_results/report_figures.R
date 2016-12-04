# In this script I reuse customized snippets from other scripts to build figures for the final report
# Cyril Matthey-Doret
# Thu Dec  1 19:14:34 2016 ------------------------------

library(knitr);library(ggplot2)

########################
# Figure 1: Expression #
########################
#path <- "/Users/cmatthe5/Documents/First_step/data/"
path <- "/home/cyril/Documents/Master/sem_1/First_step/data/"
#path <- "/home/cyril/Documents/First_step/data/"
setwd(path)
cell_lines <- c("GM12878")  # Cell lines to plot
test_lines <- c("GM12878")  # Cell lines to include in the table with corr statistics
whole_expr <- read.table("expression/enhancer_promoter/whole_exp.txt",header=T)
colnames(whole_expr) <- c("gene_id","GM12878","HUVEC","K562","NHEK","promoter","enhancer","gentype")

# Loading lincRNAs sets
e.np_linc <- read.table("enhancer_bound/all_combinations/e.np_linc_prb.bed")  # Overlap enhancer marks, but no promoter marks
ne.np_linc <- read.table("enhancer_bound/all_combinations/ne.np_linc_prb.bed")  # Overlap neither promoter marks, nor enhancer marks

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
  if(i<0.1){p <- "."}
  if(i<0.05){p <- "*"}
  if(i<0.01){p <- "**"}
  if(i<0.001){p <- "***"}
  else{p <- "-"}
  starcode <-append(p,starcode)}
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
ggplot(comp_lines)+geom_boxplot(aes(x=paste(cell.line,gentype,sep="_"),y=log10(expression),fill=cell.line), notch=T)+
  scale_x_discrete(labels=rep(gennames,length(test_lines))) +
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


# Figure 2: sequence conservation

library(ggplot2);library(gridExtra);library(plyr)

whole_cons <- read.table("../data/seq_conserv/enhancer_promoter/whole_cons.txt", header=T)
whole_cons <- whole_cons[whole_cons$promoter=="-",]
ez_class <- c()
rownames(whole_cons) <- NULL
for(r in rownames(whole_cons)){
  r <- as.numeric(r)
  if(whole_cons$gentype[r]=="pc"){
    ez_class <- append(ez_class,"PCG")
  }else{
      ez_class <- append(ez_class,ifelse(whole_cons$enhancer[r]=="+","elincRNA","other lincRNA"))
    }
}
whole_cons$gentype <- ez_class
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
  if(i<0.05){p <- "*"}
  if(i<0.01){p <- "**"}
  if(i<0.001){p <- "***"}
  else{p <- "-"}
  starcode <-append(p,starcode)}
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
ggplot(data=whole_cons)+
  facet_grid(~gr,labeller = labeller(gr = as_labeller(group_names)))+
  geom_boxplot(aes(x=gentype,y=avg_score,fill=gentype),notch=T)+
  geom_hline(data= AR_medcons,aes(yintercept = AR_medcons$val),size=1,col="#66bb66",show.legend= F)+
  scale_fill_manual(values = c("#bbbb66","#bb6666","#6666bb"))+
  #scale_linetype_manual("Title", values = 2) +
  #guides(fill = guide_legend(title = element_blank(),override.aes=list(linetype=0)))+
  #geom_text(data=wilcox_p, aes(x=c(1.5, 2), y=1, label=paste0("p=",pval)), 
  #          colour="black", inherit.aes=T, parse=T)+
  geom_label(data=med.fac, aes(x=gentype, y=x+0.01, label=x), 
             colour="black", inherit.aes=FALSE, parse=FALSE,size=3)+
  theme_bw()+ ylab("averaged phastCons")+xlab("")+guides(fill=F,line=F)+coord_cartesian(ylim=c(0,1.1))+theme(legend.position="none")+
  geom_text(data=wilcox_p,aes(x=rep(c(0.5,1),length(test_lines)*2)+increm, y=rep(c(0.94,1.01),length(test_lines)*2),label=starcode,size=5))+
  geom_line(data = df1[1:5,], aes(x = a, y = b)) +
  geom_line(data = df2[1:4,], aes(x = a, y = b)) +
  geom_line(data = df1[6:10,], aes(x = a, y = b)) +
  geom_line(data = df2[5:8,], aes(x = a, y = b)) +
  geom_line(data = df1[11:15,], aes(x = a, y = b)) +
  geom_line(data = df2[9:12,], aes(x = a, y = b)) +
  geom_line(data = df1[16:20,], aes(x = a, y = b)) +
  geom_line(data = df2[13:16,], aes(x = a, y = b))

#grid.arrange(layout_matrix=matrix(c(1,4,4,2,4,4,3,4,4),nrow = 3,byrow = T),grobs = list(dl,dp,dr,c))

# Figure 3+4: Enrichment at TAD boundaries and loop anchor
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


ggplot(data=whole_hicbound,aes(x=nicenames[c(1,3,2,4)],y=fold,fill=segment))+
  geom_bar(stat = 'identity')+
  #scale_fill_continuous(high="#DD5555",low = "#55DD55")+
  theme_bw()+guides(fill=F)+ylab("Fold enrichment")+
  #geom_hline(data=whole_hicbound[1:2,],aes(yintercept=fold[c(1,2)]))+
  xlab("")+ggtitle("Enrichment at TAD boundaries")+
  geom_text(aes(x=c(1,3,2,4),y=fold-0.1,label=fold))+
  geom_text(aes(x=c(1,3,2,4),y=fold-0.2,label=paste0("q=",qval)))+
  scale_x_discrete(labels=nicenames[c(1,3,2,4)])

ggplot(data=whole_anchors,aes(x=nicenames[c(1,3,2,4)],y=fold,fill=segment))+
  geom_bar(stat = 'identity')+
  #scale_fill_continuous(high="#DD5555",low = "#55DD55")+
  theme_bw()+ylab("Fold enrichment")+guides(fill=F)+
  #geom_hline(data=whole_anchors[1:2,],aes(yintercept=fold[c(1,2)]))+
  xlab("")+ggtitle("Enrichment at loops anchors")+
  geom_text(aes(x=c(1,3,2,4),y=fold-0.1,label=fold))+
  geom_text(aes(x=c(1,3,2,4),y=fold-0.2,label=paste0("q=",qval)))+
  scale_x_discrete(labels=nicenames[c(1,3,2,4)])

whole_bins <- read.table("GAT/out/whole_seg_10kgat_results.txt",header=T)



# Figure 4: Enrichment in architect. prots.

whole_exclu <- read.table("GAT/out/whole_seg_10kgat_exclu_insul_bs.txt",header=T)
whole_exclu <- whole_exclu[whole_exclu$track!="merged",]
whole_exclu <- whole_exclu[whole_exclu$annotation %in% c("e.np", "ne.np"),]

ggplot(data=whole_exclu[whole_exclu$segment=="exclu.CTCF",],aes(x=nicenames[c(1,3,2,4)],y=fold,fill=annotation))+
  geom_bar(stat = 'identity')+
  #scale_fill_continuous(high="#DD5555",low = "#55DD55")+
  theme_bw()+ylab("Fold enrichment")+guides(fill=F)+
  #geom_hline(data=whole_anchors[1:2,],aes(yintercept=fold[c(1,2)]))+
  xlab("")+ggtitle("Enrichment of CTCF exclusive binding sites")+
  geom_text(aes(x=c(1,3,2,4),y=fold-0.1,label=fold))+
  geom_text(aes(x=c(1,3,2,4),y=fold-0.2,label=paste0("q=",qval)))+
  scale_x_discrete(labels=nicenames[c(1,3,2,4)])

ggplot(data=whole_exclu[whole_exclu$segment=="exclu.cohesin",],aes(x=nicenames[c(1,3,2,4)],y=fold,fill=annotation))+
  geom_bar(stat = 'identity')+
  #scale_fill_continuous(high="#DD5555",low = "#55DD55")+
  theme_bw()+ylab("Fold enrichment")+guides(fill=F)+
  #geom_hline(data=whole_anchors[1:2,],aes(yintercept=fold[c(1,2)]))+
  xlab("")+ggtitle("Enrichment of cohesin exclusive binding sites")+
  geom_text(aes(x=c(1,3,2,4),y=fold-0.1,label=fold))+
  geom_text(aes(x=c(1,3,2,4),y=fold-0.2,label=paste0("q=",qval)))+
  scale_x_discrete(labels=nicenames[c(1,3,2,4)])


# Figure 6: High contact
