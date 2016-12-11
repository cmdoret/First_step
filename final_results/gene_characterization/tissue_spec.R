# This script intends to compute tissue specificity for all lincRNAs and protein-coding and store the output in a summarized table.
# Tissue specificity index (Tau) is used as a metric. I follow the methods from Kryuchkova & Robinson-Rechavi (2015).
# Cyril Matthey-Doret
# 18.11.2016

#================================================
#Loading data
options(scipen=999)  
setwd("/home/cyril/Documents/First_step/data/")
id_convert <- read.table("tissue_specificity/xloc2ensg")  # Used to convert genes ID from XLOC to ENS
colnames(id_convert) <- c("XLOC_ID","ENS_ID")
linc_RNA_tissue <- read.table("tissue_specificity/all.lincRNA.tissue.encode.txt", header=T,stringsAsFactors = F)
linc_RNA_tissue <- na.omit(linc_RNA_tissue)
linc_RNA_tissue[,2:637] <- sapply(linc_RNA_tissue[,2:637],as.numeric)
pc_tissue <- read.table("tissue_specificity/all.pcgene.tissue.encode.txt",header=T,stringsAsFactors = F)

#convert IDs to ENSG:
linc_RNA_tissue <- merge(x=linc_RNA_tissue, y= id_convert,by.x="ID",by.y="XLOC_ID",all=F)
linc_RNA_tissue <- linc_RNA_tissue[,-1]
linc_RNA_tissue <- linc_RNA_tissue[,c(637,1:636)]
pc_tissue <- merge(x=pc_tissue, y=id_convert, by.x="ID", by.y="XLOC_ID", all=F)  # merging dataframes to keep only genes that have both IDs
pc_tissue <- pc_tissue[,-1] # Removing columns containing old IDs
pc_tissue <- pc_tissue[,c(637,1:636)]  # putting column with new ID first

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

tau_pcgene <-data.frame(ID=names(taus.pc_tissue),tau=unname(taus.pc_tissue))
tau_lincRNA <-data.frame(ID=names(taus.linc_RNA_tissue),tau=unname(taus.linc_RNA_tissue))
par(mfrow=c(1,2))
hist(tau_lincRNA$tau);hist(tau_pcgene$tau)

#========================================================
# Splitting into TAD-bound and nonTAD-bound

# Loading lincRNAs sets
ne_linc <- read.table("enhancer_bound/all_combinations/ne_linc_pr.bed")  # Overlap no enhancer marks, does not take promoters into account
e_linc <- read.table("enhancer_bound/all_combinations/e_linc_pr.bed")  # Overlap enhancer marks, does not take promoters into account
ne.p_linc <- read.table("enhancer_bound/all_combinations/ne.p_linc_pr.bed")  # Overlap promoter marks, but no enhancer marks
e.np_linc <- read.table("enhancer_bound/all_combinations/e.np_linc_pr.bed")  # Overlap enhancer marks, but no promoter marks
e.p_linc <- read.table("enhancer_bound/all_combinations/e.p_linc_pr.bed")  # Overlap promoter marks, and enhancer marks
ne.np_linc <- read.table("enhancer_bound/all_combinations/ne.np_linc_pr.bed")  # Overlap neither promoter marks, nor enhancer marks

# Loading pcgenes sets
ne_pc <- read.table("enhancer_bound/all_combinations/ne_pc_pr.bed")  # Overlap no enhancer marks, does not take promoters into account
e_pc <- read.table("enhancer_bound/all_combinations/e_pc_pr.bed")  # Overlap enhancer marks, does not take promoters into account
ne.p_pc <- read.table("enhancer_bound/all_combinations/ne.p_pc_pr.bed")  # Overlap promoter marks, but no enhancer marks
e.np_pc <- read.table("enhancer_bound/all_combinations/e.np_pc_pr.bed")  # Overlap enhancer marks, but no promoter marks
e.p_pc <- read.table("enhancer_bound/all_combinations/e.p_pc_pr.bed")  # Overlap promoter marks, and enhancer marks
ne.np_pc <- read.table("enhancer_bound/all_combinations/ne.np_pc_pr.bed")  # Overlap neither promoter marks, nor enhancer marks

#adding colnames
colnames(ne_linc)= colnames(e_linc)=colnames(ne.p_linc)=colnames(e.np_linc)=colnames(e.p_linc)=
  colnames(ne.np_linc)=colnames(ne_pc)=colnames(e_pc)=colnames(ne.p_pc)=colnames(e.np_pc)=colnames(e.p_pc)=
  colnames(ne.np_pc) <-c("chr", "start", "end", "gene", "strand")


#splitting specificity values into categories

crop <- function(s){sub("\\..*","",s)}
tau_ne_linc <-tau_lincRNA[tau_lincRNA$ID %in% crop(ne_linc$gene),]
tau_e_linc <-tau_lincRNA[tau_lincRNA$ID %in% crop(e_linc$gene),]
tau_e.p_linc <-tau_lincRNA[tau_lincRNA$ID %in% crop(e.p_linc$gene),]
tau_ne.p_linc <-tau_lincRNA[tau_lincRNA$ID %in% crop(ne.p_linc$gene),]
tau_e.np_linc <-tau_lincRNA[tau_lincRNA$ID %in% crop(e.np_linc$gene),]
tau_ne.np_linc <-tau_lincRNA[tau_lincRNA$ID %in% crop(ne.np_linc$gene),]

tau_ne_pc <-tau_pcgene[tau_pcgene$ID %in% crop(ne_pc$gene),]
tau_e_pc <-tau_pcgene[tau_pcgene$ID %in% crop(e_pc$gene),]
tau_e.p_pc <-tau_pcgene[tau_pcgene$ID %in% crop(e.p_pc$gene),]
tau_ne.p_pc <-tau_pcgene[tau_pcgene$ID %in% crop(ne.p_pc$gene),]
tau_e.np_pc <-tau_pcgene[tau_pcgene$ID %in% crop(e.np_pc$gene),]
tau_ne.np_pc <-tau_pcgene[tau_pcgene$ID %in% crop(ne.np_pc$gene),]

# Putting all values in the same dataframe, for the sake of convenience
whole_tau <- rbind(cbind(tau_e_linc,promoter=rep("both"),enhancer=rep("+"),gentype=rep("linc")), 
                   cbind(tau_ne_linc,promoter=rep("both"),enhancer=rep("-"),gentype=rep("linc")),
                   cbind(tau_e.np_linc,promoter=rep("-"),enhancer=rep("+"),gentype=rep("linc")),
                   cbind(tau_ne.p_linc,promoter=rep("+"),enhancer=rep("-"),gentype=rep("linc")),
                   cbind(tau_e.p_linc,promoter=rep("+"),enhancer=rep("+"),gentype=rep("linc")),
                   cbind(tau_ne.np_linc,promoter=rep("-"),enhancer=rep("-"),gentype=rep("linc")),
                   cbind(tau_e_pc,promoter=rep("both"),enhancer=rep("+"),gentype=rep("pc")), 
                   cbind(tau_ne_pc,promoter=rep("both"),enhancer=rep("-"),gentype=rep("pc")),
                   cbind(tau_e.np_pc,promoter=rep("-"),enhancer=rep("+"),gentype=rep("pc")),
                   cbind(tau_ne.p_pc,promoter=rep("+"),enhancer=rep("-"),gentype=rep("pc")),
                   cbind(tau_e.p_pc,promoter=rep("+"),enhancer=rep("+"),gentype=rep("pc")),
                   cbind(tau_ne.np_pc,promoter=rep("-"),enhancer=rep("-"),gentype=rep("pc")))
write.table(x = whole_tau,file = "tissue_specificity/enhancer_promoter/whole_tau.txt",quote = F,sep = "\t",row.names = F,col.names = T)
