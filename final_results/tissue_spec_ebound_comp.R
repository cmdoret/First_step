# This script intends to compare tissue specificity in expression (in RPKM) between ebound and pbound lincRNAs and protein coding genes.
# Cyril Matthey-DOret
# 07.10.2016

#================================================
#Loading data
options(scipen=999)
setwd("/Users/cmatthe5/Documents/First_step/data/")
#setwd("/home/cyril/Documents/First_step/data/")
id_convert <- read.table("tissue_specificity/xloc2ensg")
colnames(id_convert) <- c("XLOC_ID","ENS_ID")
linc_RNA_tissue <- read.table("tissue_specificity/all.lincRNA.tissue.encode.txt", header=T,stringsAsFactors = F)
linc_RNA_tissue <- na.omit(linc_RNA_tissue)
linc_RNA_tissue[,2:637] <- sapply(linc_RNA_tissue[,2:637],as.numeric)
pc_tissue <- read.table("tissue_specificity/all.pcgene.tissue.encode.txt",header=T,stringsAsFactors = F)

#convert IDs to ENSG:
linc_RNA_tissue <- merge(x=linc_RNA_tissue, y= id_convert,by.x="ID",by.y="XLOC_ID",all=F)
linc_RNA_tissue <- linc_RNA_tissue[,-1]
linc_RNA_tissue <- linc_RNA_tissue[,c(637,1:636)]
pc_tissue <- merge(x=pc_tissue, y=id_convert, by.x="ID", by.y="XLOC_ID", all=F)
pc_tissue <- pc_tissue[,-1]
pc_tissue <- pc_tissue[,c(637,1:636)]

#Calculation of tissue specificity:
### following methods from Robinson's TS review ###
## http://bib.oxfordjournals.org/content/early/2016/02/17/bib.bbw008.full


calcTDR <- function( mat, cutoff=0 ) {  
  # transcript detection rate: in what proportion of samples were each transcript seen?
  # cutoff optional
  100 * rowSums(mat>cutoff) / ncol(mat)
}
calcTau <- function( v, cutoff=0 ) {     
  # Tau measure of specificity (0: generic, 1:specific)
  # cutoff optional
  v[v<cutoff] <- cutoff
  sum(1-(v/max(v)))/(length(v)-1)
}
calcNormExpr <- function(mat, cutoff=0.1) { 
  # transform: shifted log10
  # cutoff optional
  mat[mat<=cutoff] <- 0
  log10( mat + 1 )
}

# transform expression

rownames(linc_RNA_tissue)<- linc_RNA_tissue$ENS_ID
linc_RNA_tissue.norm <- calcNormExpr(linc_RNA_tissue[,-1])
rownames(pc_tissue)<- pc_tissue$ENS_ID
pc_tissue.norm <- calcNormExpr(pc_tissue[,-1])

# get "naive" TDRs to filter transcripts never detected
pc_tdr <- calcTDR(pc_tissue.norm) 
linc_RNA_tdr <- calcTDR(linc_RNA_tissue.norm) 

linc_RNA_tissue.norm_filt <- linc_RNA_tissue.norm[linc_RNA_tdr>0,]
pc_tissue.norm_filt <- pc_tissue.norm[pc_tdr>0,]


# calculate tissue-averaged Taus
taus.pc_tissue <- apply(pc_tissue.norm_filt, MARGIN=1, FUN=calcTau)
taus.linc_RNA_tissue <- apply(linc_RNA_tissue.norm_filt, MARGIN=1, FUN=calcTau)

tau_pc <-data.frame(ID=names(taus.pc_tissue),tau=unname(taus.pc_tissue))
tau_linc_RNA <-data.frame(ID=names(taus.linc_RNA_tissue),tau=unname(taus.linc_RNA_tissue))
par(mfrow=c(1,2))
hist(tau_linc_RNA$tau);hist(tau_pc$tau)
# lincRNA are very tissue-specific.
#========================================================
# Splitting into prom.assoc. and enh.assoc.

#loading genes sets
elinc <- read.table("enhancer_bound/elinc_prb.bed")
plinc<- read.table("enhancer_bound/plinc_prb.bed")
epc<- read.table("enhancer_bound/epc_prb.bed")
ppc<- read.table("enhancer_bound/ppc_prb.bed")

# adding colnames
colnames(elinc)= colnames(plinc)=colnames(epc)=colnames(ppc) <-c("chr", "start", "end", "gene", "strand")

#splitting specificity values into TADb and nTADb
crop <- function(s){sub("\\..*","",s)} # removing version in gene names

tau_elinc <-tau_linc_RNA[tau_linc_RNA$ID %in% crop(elinc$gene),]
tau_plinc <-tau_linc_RNA[tau_linc_RNA$ID %in% crop(plinc$gene),]
tau_epc <-tau_pc[tau_pc$ID %in% crop(epc$gene),]
tau_ppc <-tau_pc[tau_pc$ID %in% crop(ppc$gene),]


# Putting all values in the same dataframe, for the sake of convenience
whole_tau <- rbind(cbind(tau_elinc,assoc=rep("e"),gentype=rep("lincRNA")), 
                   cbind(tau_plinc,assoc=rep("p"),gentype=rep("lincRNA")), 
                   cbind(tau_epc,assoc=rep("e"),gentype=rep("pc")), 
                   cbind(tau_ppc,assoc=rep("p"),gentype=rep("pc")))
write.table(x = whole_tau,file = "tissue_specificity/enhancer_bound/whole_tau.txt",quote = F,sep = "\t",row.names = F,col.names = T)
#========================

#Visualizing:

library(plyr)
short_med<-function(x){return(round(median(log10(x)),3))}
short_wilcox <- function(x,y){return(format(wilcox.test(x, y)$p.value,digits=3))}

#dataf frame for annotation of median
med.fac = ddply(whole_tau, .(gentype, assoc), function(.d)
  data.frame(x=median(round(.d$tau,3))))
#data frame for annotation of p-value
pl<-c()
gl<-c()


for(g in levels(whole_tau$gentype)){
  p <- short_wilcox(whole_tau[whole_tau$gentype==g &  whole_tau$assoc=="e","tau"],
                    whole_tau[whole_tau$gentype==g & whole_tau$assoc=="p","tau"])
  pl<-append(pl,p)
  gl<-append(gl,g)
}

wilcox_p <- data.frame(pval=pl,gentype=gl)
library(ggplot2)
ggplot(data=whole_tau)+
  facet_grid(.~gentype)+
  geom_boxplot(aes(x=assoc, y=tau,fill=assoc),notch=T)+
  geom_text(data=med.fac, aes(x=assoc, y=0, label=x), 
            colour="black", inherit.aes=FALSE, parse=FALSE)+
  geom_text(data=wilcox_p, aes(x=1.5, y=0.25, label=paste0("p-value=",pval)), 
            colour="black", inherit.aes=FALSE, parse=FALSE)+
  theme_bw()
