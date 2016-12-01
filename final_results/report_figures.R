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


increm <- c()
for(i in seq(1,3*length(test_lines),3)){increm <- append(increm,rep(i,2))}  # building vector for positioning p-values between boxes
ggplot(comp_lines)+geom_boxplot(aes(x=paste(cell.line,gentype,sep="_"),y=log10(expression),fill=cell.line))+
  scale_x_discrete(labels=rep(gennames,length(test_lines))) +
  theme_bw()+xlab("Gene class")+ylab("Log10 expression (FPKM)")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  annotate(geom = 'text',x=rep(c(0.5,1),length(test_lines))+increm, y=rep(c(5.7,6.7),length(test_lines)),label=wilcox_p$starcode,size=5)+
  geom_line(data = df1[1:5,], aes(x = a, y = b)) +
  geom_line(data = df2[1:4,], aes(x = a, y = b)) +
  geom_line(data = df1[6:10,], aes(x = a, y = b)) +
  geom_line(data = df2[5:8,], aes(x = a, y = b)) +
  geom_line(data = df1[11:15,], aes(x = a, y = b)) +
  geom_line(data = df2[9:12,], aes(x = a, y = b)) +
  geom_line(data = df1[16:20,], aes(x = a, y = b)) +
  geom_line(data = df2[13:16,], aes(x = a, y = b))
# That's a mess, but it works
