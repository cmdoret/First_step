# This script intends to compare tissue specificity in expression (in RPKM) between TADbound and non-TADbound lincRNAs and protein coding genes.
# Cyril Matthey-DOret
# 07.10.2016

#================================================
#Loading data
options(scipen=999)
setwd("/home/cyril/Documents/First_step/data/")
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
# Splitting into TAD-bound and nonTAD-bound

#loading TAD-bound lincRNAs sets
Tb_lincRNA5 <- read.table("linc_RNA/merged/TADbound-lincRNA5.bed")
Tb_lincRNA10 <- read.table("linc_RNA/merged/TADbound-lincRNA10.bed")
Tb_lincRNA20 <- read.table("linc_RNA/merged/TADbound-lincRNA20.bed")
#loading TAD-bound pcgenes sets
Tb_pc5 <- read.table("pc_genes/merged/TADbound-pcgene5.bed")
Tb_pc10 <- read.table("pc_genes/merged/TADbound-pcgene10.bed")
Tb_pc20 <- read.table("pc_genes/merged/TADbound-pcgene20.bed")
#loading non-TAD-bound lincRNAs sets
nTb_lincRNA5 <- read.table("linc_RNA/merged/nonTADbound-lincRNA5.bed")
nTb_lincRNA10 <- read.table("linc_RNA/merged/nonTADbound-lincRNA10.bed")
nTb_lincRNA20 <- read.table("linc_RNA/merged/nonTADbound-lincRNA20.bed")
#loading non-TAD-bound pcgenes sets
nTb_pc5 <- read.table("pc_genes/merged/nonTADbound-pcgene5.bed")
nTb_pc10 <- read.table("pc_genes/merged/nonTADbound-pcgene10.bed")
nTb_pc20 <- read.table("pc_genes/merged/nonTADbound-pcgene20.bed")

# adding colnames
colnames(Tb_lincRNA5)= colnames(Tb_lincRNA10)=colnames(Tb_lincRNA20)=colnames(nTb_lincRNA5)=colnames(nTb_lincRNA10)=
  colnames(nTb_lincRNA20)=colnames(Tb_pc5)=colnames(Tb_pc10)=colnames(Tb_pc20)=colnames(nTb_pc5)=colnames(nTb_pc10)=
  colnames(nTb_pc20) <-c("chr", "start", "end", "gene", "strand")

#splitting specificity values into TADb and nTADb
crop <- function(s){sub("\\..*","",s)}
tau_Tb_lincRNA5 <-tau_linc_RNA[tau_linc_RNA$ID %in% crop(Tb_lincRNA5$gene),]
tau_Tb_lincRNA10 <-tau_linc_RNA[tau_linc_RNA$ID %in% crop(Tb_lincRNA10$gene),]
tau_Tb_lincRNA20 <-tau_linc_RNA[tau_linc_RNA$ID %in% crop(Tb_lincRNA20$gene),]

tau_nTb_lincRNA5 <-tau_linc_RNA[tau_linc_RNA$ID %in% crop(nTb_lincRNA5$gene),]
tau_nTb_lincRNA10 <-tau_linc_RNA[tau_linc_RNA$ID %in% crop(nTb_lincRNA10$gene),]
tau_nTb_lincRNA20 <-tau_linc_RNA[tau_linc_RNA$ID %in% crop(nTb_lincRNA20$gene),]

tau_Tb_pc5 <-tau_pc[tau_pc$ID %in% crop(Tb_pc5$gene),]
tau_Tb_pc10 <-tau_pc[tau_pc$ID %in% crop(Tb_pc10$gene),]
tau_Tb_pc20 <-tau_pc[tau_pc$ID %in% crop(Tb_pc20$gene),]

tau_nTb_pc5 <-tau_pc[tau_pc$ID %in% crop(nTb_pc5$gene),]
tau_nTb_pc10 <-tau_pc[tau_pc$ID %in% crop(nTb_pc10$gene),]
tau_nTb_pc20 <-tau_pc[tau_pc$ID %in% crop(nTb_pc20$gene),]

#writing into bed files for further use:
write.table(tau_Tb_lincRNA5,file = "tissue_specificity/merged/tau_Tb_lincRNA5.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(tau_Tb_lincRNA10,file = "tissue_specificity/merged/tau_Tb_lincRNA10.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(tau_Tb_lincRNA20,file = "tissue_specificity/merged/tau_Tb_lincRNA20.txt",sep="\t",quote = F,col.names = F,row.names = F)

write.table(tau_nTb_lincRNA5,file = "tissue_specificity/merged/tau_nTb_lincRNA5.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(tau_nTb_lincRNA10,file = "tissue_specificity/merged/tau_nTb_lincRNA10.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(tau_nTb_lincRNA20,file = "tissue_specificity/merged/tau_nTb_lincRNA20.txt",sep="\t",quote = F,col.names = F,row.names = F)

write.table(tau_Tb_pc5,file = "tissue_specificity/merged/tau_Tb_pc5.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(tau_Tb_pc10,file = "tissue_specificity/merged/tau_Tb_pc10.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(tau_Tb_pc20,file = "tissue_specificity/merged/tau_Tb_pc20.txt",sep="\t",quote = F,col.names = F,row.names = F)

write.table(tau_nTb_pc5,file = "tissue_specificity/merged/tau_nTb_pc5.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(tau_nTb_pc10,file = "tissue_specificity/merged/tau_nTb_pc10.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(tau_nTb_pc20,file = "tissue_specificity/merged/tau_nTb_pc20.txt",sep="\t",quote = F,col.names = F,row.names = F)

# Putting all values in the same dataframe, for the sake of convenience
whole_tau <- rbind(cbind(tau_Tb_pc5,threshold=rep("5"),TAD=rep("Tb"),gentype=rep("pc")), 
                   cbind(tau_Tb_pc10,threshold=rep("10"),TAD=rep("Tb"),gentype=rep("pc")), 
                   cbind(tau_Tb_pc20,threshold=rep("20"),TAD=rep("Tb"),gentype=rep("pc")), 
                   cbind(tau_nTb_pc5,threshold=rep("5"),TAD=rep("nTb"),gentype=rep("pc")), 
                   cbind(tau_nTb_pc10,threshold=rep("10"),TAD=rep("nTb"),gentype=rep("pc")), 
                   cbind(tau_nTb_pc20,threshold=rep("20"),TAD=rep("nTb"),gentype=rep("pc")),
                   cbind(tau_Tb_lincRNA5,threshold=rep("5"),TAD=rep("Tb"),gentype=rep("lincRNA")), 
                   cbind(tau_Tb_lincRNA10,threshold=rep("10"),TAD=rep("Tb"),gentype=rep("lincRNA")), 
                   cbind(tau_Tb_lincRNA20,threshold=rep("20"),TAD=rep("Tb"),gentype=rep("lincRNA")), 
                   cbind(tau_nTb_lincRNA5,threshold=rep("5"),TAD=rep("nTb"),gentype=rep("lincRNA")), 
                   cbind(tau_nTb_lincRNA10,threshold=rep("10"),TAD=rep("nTb"),gentype=rep("lincRNA")), 
                   cbind(tau_nTb_lincRNA20,threshold=rep("20"),TAD=rep("nTb"),gentype=rep("lincRNA")))
write.table(x = whole_tau,file = "tissue_specificity/merged/whole_tau.txt",quote = F,sep = "\t",row.names = F,col.names = T)
#Visualizing:

library(plyr)
short_med<-function(x){return(round(median(log10(x)),3))}
short_wilcox <- function(x,y){return(format(wilcox.test(x, y)$p.value,digits=3))}

#dataf frame for annotation of median
med.fac = ddply(whole_tau, .(gentype, threshold,TAD), function(.d)
  data.frame(x=median(round(.d$tau,3))))
#data frame for annotation of p-value
pl<-c()
gl<-c()
tl<-c()
for(t in levels(whole_tau$threshold)){
  for(g in levels(whole_tau$gentype)){
    p <- short_wilcox(whole_tau[whole_tau$threshold==t & whole_tau$gentype==g &  whole_tau$TAD=="Tb","tau"],
                      whole_tau[whole_tau$threshold==t &  whole_tau$gentype==g & whole_tau$TAD=="nTb","tau"])
    pl<-append(pl,p)
    tl<-append(tl,t)
    gl<-append(gl,g)
  }
}
wilcox_p <- data.frame(pval=pl,threshold=tl,gentype=gl)
library(ggplot2)
ggplot(data=whole_tau)+
  facet_grid(gentype~threshold)+
  geom_boxplot(aes(x=TAD, y=tau,fill=TAD),notch=T)+
  geom_text(data=med.fac, aes(x=TAD, y=0, label=x), 
            colour="black", inherit.aes=FALSE, parse=FALSE)+
  geom_text(data=wilcox_p, aes(x=1.5, y=0.25, label=paste0("p-value=",pval)), 
            colour="black", inherit.aes=FALSE, parse=FALSE)+
  theme_bw()
