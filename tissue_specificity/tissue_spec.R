# This script intends to compare tissue specificity in expression (in RPKM) between TADbound and non-TADbound lincRNAs and protein coding genes.
# Cyril Matthey-DOret
# 07.10.2016

#================================================
#Loading data
options(scipen=999)
setwd("/home/cyril/Documents/First_step/data/tissue_specificity/")
id_convert <- read.table("xloc2ensg")
colnames(id_convert) <- c("XLOC_ID","ENS_ID")
linc_RNA_tissue <- read.table("all.lincRNA.tissue.encode.txt", header=T,stringsAsFactors = F)
linc_RNA_tissue <- na.omit(linc_RNA_tissue)
linc_RNA_tissue[,2:637] <- sapply(linc_RNA_tissue[,2:637],as.numeric)
pc_tissue <- read.table("all.pcgene.tissue.encode.txt",header=T,stringsAsFactors = F)

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
linc_RNA_tissue.norm <- calcNormExpr(linc_RNA_tissue[,-1])
pc_tissue.norm <- calcNormExpr(pc_tissue[,-1])

# get "naive" TDRs to filter transcripts never detected
pc_tdr <- calcTDR(pc_tissue.norm) 
linc_RNA_tdr <- calcTDR(linc_RNA_tissue.norm) 

linc_RNA_tissue.norm_filt <- linc_RNA_tissue.norm[linc_RNA_tdr>0,]
pc_tissue.norm_filt <- pc_tissue.norm[pc_tdr>0,]


# calculate tissue-averaged Taus
taus.pc_tissue <- apply(pc_tissue.norm_filt, MARGIN=1, FUN=calcTau)
taus.linc_RNA_tissue <- apply(linc_RNA_tissue.norm_filt, MARGIN=1, FUN=calcTau)

tau_pc <-merge(x=pc_tissue$ENS_ID,y=unname(taus.pc_tissue),by.x=rownames(pc_tissue),by.y=names(taus.pc_tissue))
hist(unname(taus.pc_tissue))
