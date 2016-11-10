# In this script I check if there are correlations between the expression of different gene types and the amount of 
# DNA-DNA contacts at their location. This analysis is done across 4 different cell types.
# Cyril Matthey-Doret
# Thu Nov 10 13:58:33 2016 ------------------------------

setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")

test_contact <- read.table("TAD_contact/all.lincRNA.GM12878.HiC.contact.txt", header=T)
test_expr <- read.table("expression/LCL.lincRNA.expression.txt")
colnames(test_expr) <- c("gene", "expression")


contexpr <- function(hic, expr, cell_line="", gtype=""){  # This function computes the correlation between expression and DNA contacts for a set of genes.
# Both arguments need to be a list of genes containing the gene ID and its median expression level/amount of DNA-DNA contact.
  merged <- merge(x = hic, y = expr, by = "gene", all = F)
  results <- c(cor.test(merged$diam, merged$expression, method="spearman")$estimate, 
               cor.test(merged$diam, merged$expression, method="spearman")$p.value)
  smoothScatter(log10(merged$diam), log10(merged$expression), 
                xlab = "DNA-DNA interactions",ylab="Median expression levels", main=paste(cell_line,gtype,sep=": "))
  text(x = c(6.5,6.5),y=c(1.5,1), labels=c(paste0("rho = ", round(results[1],3)), paste0("p-value = ", round(results[2],3))))
  
}



#================================
# Analysis for elincRNAs and nelincRNAs

elinc <- read.table("enhancer_bound/elinc_prb.bed")
nelinc <- read.table("enhancer_bound/nelinc_prb.bed")
colnames(elinc) = colnames(nelinc) <- c("chr", "start", "end", "gene", "strand")
par(mfrow=c(2,2))
#for(c in c("GM12878", "HUVEC", "K562", "NHEK")){
for(c in c("GM12878")){
  linc_contact <- read.table(paste0("TAD_contact/all.lincRNA.",c,".HiC.contact.txt"), header=T)
  linc_expr <- read.table("expression/LCL.lincRNA.expression.txt")
  colnames(linc_expr) <- c("gene", "expression")
  elinc_expr <- linc_expr[linc_expr$gene %in% elinc$gene,]
  nelinc_expr <-linc_expr[linc_expr$gene %in% nelinc$gene,]
  elinc_contact <- linc_contact[linc_contact$gene %in% elinc$gene,]
  nelinc_contact <- linc_contact[linc_contact$gene %in% nelinc$gene,]
  contexpr(elinc_contact, elinc_expr,cell_line=c,gtype="elincRNA")
  contexpr(nelinc_contact, nelinc_expr,cell_line=c,gtype="nelincRNA")
}

#================================
# Analysis for epc and nepc

epc <- read.table("enhancer_bound/epc_prb.bed")
nepc <- read.table("enhancer_bound/nepc_prb.bed")
colnames(epc) = colnames(nepc) <- c("chr", "start", "end", "gene", "strand")
#par(mfrow=c(2,2))
#for(c in c("GM12878", "HUVEC", "K562", "NHEK")){
for(c in c("GM12878")){
  pc_contact <- read.table(paste0("TAD_contact/all.pcgene.",c,".HiC.contact.txt"), header=T)
  pc_expr <- read.table("expression/LCL.pcgene.expression.txt")
  colnames(pc_expr) <- c("gene", "expression")
  epc_expr <- pc_expr[pc_expr$gene %in% epc$gene,]
  nepc_expr <-pc_expr[pc_expr$gene %in% nepc$gene,]
  epc_contact <- pc_contact[pc_contact$gene %in% epc$gene,]
  nepc_contact <- pc_contact[pc_contact$gene %in% nepc$gene,]
  contexpr(epc_contact, epc_expr,cell_line=c,gtype="epc gene")
  contexpr(nepc_contact, nepc_expr,cell_line=c,gtype="nepc gene")
}

summary(glm(log10(nelinc_expr$expression[nelinc_expr$gene %in% nelinc_contact$gene]) ~ 
          nelinc_contact$diam[nelinc_contact$gene %in% nelinc_expr$gene]))
