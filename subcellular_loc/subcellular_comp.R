# The purpose of this script is to compare the sucellular localization of TADbound and non-TADbound lincRNAs.
# It also compares the subcellular localizatioin of TADbound and non-TADbound protein coding genes.

# Cyril Matthey-Doret, 05.10.2016

library(ggplot2)
library(gridExtra)

######################################################

#Loading data:
#setwd("/home/cyril/Documents/First_step/data/")
#Loading data:
setwd("/home/cyril/Documents/First_step/data/")
#setwd("/home/cyril/Documents/Master/sem_1/First step/data/")
bedcol <- c("chr", "start", "end", "gene", "strand")
loc_lincRNA <- read.table("subcell_loc/all.lincRNA.GM12878.subcellular.ratio.txt", header = T)
loc_pcgene <- read.table("subcell_loc/all.pcgene.GM12878.subcellular.ratio.txt", header = T)
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

colnames(Tb_lincRNA5)= colnames(Tb_lincRNA10)=colnames(Tb_lincRNA20)=colnames(nTb_lincRNA5)=colnames(nTb_lincRNA10)=
  colnames(nTb_lincRNA20)=colnames(Tb_pc5)=colnames(Tb_pc10)=colnames(Tb_pc20)=colnames(nTb_pc5)=colnames(nTb_pc10)=
  colnames(nTb_pc20) <-c("chr", "start", "end", "gene", "strand")


#Splitting localisation data into TADbound (Tb) and non-TADbound (nTb)

loc_Tb_lincRNA5 <-loc_lincRNA[loc_lincRNA$gene %in% Tb_lincRNA5$gene,]
loc_Tb_lincRNA10 <-loc_lincRNA[loc_lincRNA$gene %in% Tb_lincRNA10$gene,]
loc_Tb_lincRNA20 <-loc_lincRNA[loc_lincRNA$gene %in% Tb_lincRNA20$gene,]

loc_nTb_lincRNA5 <-loc_lincRNA[loc_lincRNA$gene %in% nTb_lincRNA5$gene,]
loc_nTb_lincRNA10 <-loc_lincRNA[loc_lincRNA$gene %in% nTb_lincRNA10$gene,]
loc_nTb_lincRNA20 <-loc_lincRNA[loc_lincRNA$gene %in% nTb_lincRNA20$gene,]

loc_Tb_pc5 <-loc_pcgene[loc_pcgene$gene %in% Tb_pc5$gene,]
loc_Tb_pc10 <-loc_pcgene[loc_pcgene$gene %in% Tb_pc10$gene,]
loc_Tb_pc20 <-loc_pcgene[loc_pcgene$gene %in% Tb_pc20$gene,]

loc_nTb_pc5 <-loc_pcgene[loc_pcgene$gene %in% nTb_pc5$gene,]
loc_nTb_pc10 <-loc_pcgene[loc_pcgene$gene %in% nTb_pc10$gene,]
loc_nTb_pc20 <-loc_pcgene[loc_pcgene$gene %in% nTb_pc20$gene,]

#Writing into new files
write.table(loc_Tb_lincRNA5,file = "subcell_loc/merged/loc_Tb_lincRNA5.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(loc_Tb_lincRNA10,file = "subcell_loc/merged/loc_Tb_lincRNA10.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(loc_Tb_lincRNA20,file = "subcell_loc/merged/loc_Tb_lincRNA20.txt",sep="\t",quote = F,col.names = F,row.names = F)

write.table(loc_nTb_lincRNA5,file = "subcell_loc/merged/loc_nTb_lincRNA5.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(loc_nTb_lincRNA10,file = "subcell_loc/merged/loc_nTb_lincRNA10.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(loc_nTb_lincRNA20,file = "subcell_loc/merged/loc_nTb_lincRNA20.txt",sep="\t",quote = F,col.names = F,row.names = F)

write.table(loc_Tb_pc5,file = "subcell_loc/merged/loc_Tb_pc5.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(loc_Tb_pc10,file = "subcell_loc/merged/loc_Tb_pc10.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(loc_Tb_pc20,file = "subcell_loc/merged/loc_Tb_pc20.txt",sep="\t",quote = F,col.names = F,row.names = F)

write.table(loc_nTb_pc5,file = "subcell_loc/merged/loc_nTb_pc5.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(loc_nTb_pc10,file = "subcell_loc/merged/loc_nTb_pc10.txt",sep="\t",quote = F,col.names = F,row.names = F)
write.table(loc_nTb_pc20,file = "subcell_loc/merged/loc_nTb_pc20.txt",sep="\t",quote = F,col.names = F,row.names = F)

#============================================================

#Comparing overall mean ratio between Tb and nTb:

#lincRNAs:
hist_linc <-ggplot()+
  geom_histogram(data=loc_Tb_lincRNA, aes(x=ratio, y=..density..), fill="#0000dd", alpha=0.5, bins = 300)+
  geom_histogram(data=loc_nTb_lincRNA, aes(x=ratio, y=..density..), fill="#bb0000", alpha=0.5, bins = 300)+
  ggtitle("lincRNAs")
summary(loc_Tb_lincRNA$ratio);summary(loc_nTb_lincRNA$ratio)

#Protein-coding genes:
hist_pc<-ggplot()+
  geom_histogram(data=loc_Tb_pc, aes(x=ratio, y=..density..), fill="#0000dd", alpha=0.5, bins = 300)+
  geom_histogram(data=loc_nTb_pc, aes(x=ratio, y=..density..), fill="#bb0000", alpha=0.5, bins = 300)+
  ggtitle("Protein-coding genes")
summary(loc_Tb_lincRNA$ratio);summary(loc_nTb_lincRNA$ratio)

grid.arrange(hist_linc, hist_pc)
#It doesn't seem TADbound lincRNAs/proteins are more localized in the nucleus/cytoplasm than non-TADbound ones.
boxplot(notch=T, log10(loc_Tb_lincRNA$ratio), log10(loc_nTb_lincRNA$ratio))
wilcox.test(loc_Tb_lincRNA$ratio, loc_nTb_lincRNA$ratio)
#==========================================================

#building large, single dataframe to make data more convenient.
CData <-function(df){
  return(cbind(df,fact=rep(deparse(substitute(df)),length(df[,1]))))
}
whole_loc <- rbind(CData(loc_Tb_pc5), CData(loc_Tb_pc10), CData(loc_Tb_pc20), CData(loc_nTb_pc5), CData(loc_nTb_pc10), CData(loc_nTb_pc20),
                   CData(loc_Tb_lincRNA5), CData(loc_Tb_lincRNA10), CData(loc_Tb_lincRNA20), CData(loc_nTb_lincRNA5), CData(loc_nTb_lincRNA10), CData(loc_nTb_lincRNA20))

whole_loc <- rbind(cbind(loc_Tb_pc5,threshold=rep("5"),TAD=rep("Tb"),gentype=rep("pc")), 
                   cbind(loc_Tb_pc10,threshold=rep("10"),TAD=rep("Tb"),gentype=rep("pc")), 
                   cbind(loc_Tb_pc20,threshold=rep("20"),TAD=rep("Tb"),gentype=rep("pc")), 
                   cbind(loc_nTb_pc5,threshold=rep("5"),TAD=rep("nTb"),gentype=rep("pc")), 
                   cbind(loc_nTb_pc10,threshold=rep("10"),TAD=rep("nTb"),gentype=rep("pc")), 
                   cbind(loc_nTb_pc20,threshold=rep("20"),TAD=rep("nTb"),gentype=rep("pc")),
                   cbind(loc_Tb_lincRNA5,threshold=rep("5"),TAD=rep("Tb"),gentype=rep("lincRNA")), 
                   cbind(loc_Tb_lincRNA10,threshold=rep("10"),TAD=rep("Tb"),gentype=rep("lincRNA")), 
                   cbind(loc_Tb_lincRNA20,threshold=rep("20"),TAD=rep("Tb"),gentype=rep("lincRNA")), 
                   cbind(loc_nTb_lincRNA5,threshold=rep("5"),TAD=rep("nTb"),gentype=rep("lincRNA")), 
                   cbind(loc_nTb_lincRNA10,threshold=rep("10"),TAD=rep("nTb"),gentype=rep("lincRNA")), 
                   cbind(loc_nTb_lincRNA20,threshold=rep("20"),TAD=rep("nTb"),gentype=rep("lincRNA")))
write.table(x = whole_loc,file = "subcell_loc/merged/whole_loc.txt",quote = F,sep = "\t",row.names = F,col.names = T)
library(plyr)
short_med<-function(x){return(round(median(log10(x)),3))}
short_wilcox <- function(x,y){return(format(wilcox.test(x, y)$p.value,digits=3))}

med.fac = ddply(whole_loc, .(gentype, threshold,TAD), function(.d)
  data.frame(x=round(median(log10(.d$ratio)),3)))
#data frame for annotation of p-value
pl<-c()
gl<-c()
tl<-c()
for(t in levels(whole_loc$threshold)){
  for(g in levels(whole_loc$gentype)){
    p <- short_wilcox(whole_loc[whole_loc$threshold==t & whole_loc$gentype==g &  whole_loc$TAD=="Tb","ratio"],
                      whole_loc[whole_loc$threshold==t &  whole_loc$gentype==g & whole_loc$TAD=="nTb","ratio"])
    pl<-append(pl,p)
    tl<-append(tl,t)
    gl<-append(gl,g)
  }
}
wilcox_p <- data.frame(pval=pl,threshold=tl,gentype=gl)
library(ggplot2)
ggplot(data=whole_loc)+
  facet_grid(gentype~threshold)+
  geom_boxplot(aes(x=TAD, y=ratio,fill=TAD),notch=T)+
  geom_text(data=med.fac, aes(x=TAD, y=0, label=x), 
            colour="black", inherit.aes=FALSE, parse=FALSE)+
  geom_text(data=wilcox_p, aes(x=1.5, y=0.25, label=paste0("p-value=",pval)), 
            colour="black", inherit.aes=FALSE, parse=FALSE)+
  theme_bw()


l5 <-ggplot(data=whole_loc[whole_loc$fact %in% c("loc_Tb_lincRNA5","loc_nTb_lincRNA5"),])+
  geom_boxplot(aes(fact, log10(ratio)), fill=c("darkred","darkblue"), notch = T)+theme_bw()+ggtitle("lincRNA, 5%")+
  annotate(x=c(1, 2),y=c(-1,-1),geom = "text", label=c(short_med(whole_loc$ratio[whole_loc$fact=="loc_Tb_lincRNA5"]), 
                                                       short_med(whole_loc$ratio[whole_loc$fact=="loc_nTb_lincRNA5"])))+
  scale_x_discrete(name=paste("p-value = ", short_wilcox(whole_loc$ratio[whole_loc$fact=="loc_Tb_lincRNA5"],
                                                         whole_loc$ratio[whole_loc$fact=="loc_nTb_lincRNA5"]), sep=" "),
                   labels=c("TADb", "non-TADb"))
l10 <-ggplot(data=whole_loc[whole_loc$fact %in% c("loc_Tb_lincRNA10","loc_nTb_lincRNA10"),])+
  geom_boxplot(aes(fact, log10(ratio)), fill=c("darkred","darkblue"), notch = T)+theme_bw()+ggtitle("lincRNA, 10%")+
  annotate(x=c(1, 2),y=c(-1,-1),geom = "text", label=c(short_med(whole_loc$ratio[whole_loc$fact=="loc_Tb_lincRNA10"]), 
                                                       short_med(whole_loc$ratio[whole_loc$fact=="loc_nTb_lincRNA10"])))+
  scale_x_discrete(name=paste("p-value = ", short_wilcox(whole_loc$ratio[whole_loc$fact=="loc_Tb_lincRNA10"],
                                                         whole_loc$ratio[whole_loc$fact=="loc_nTb_lincRNA10"]), sep=" "),
                   labels=c("TADb", "non-TADb"))
l20 <-ggplot(data=whole_loc[whole_loc$fact %in% c("loc_Tb_lincRNA20","loc_nTb_lincRNA20"),])+
  geom_boxplot(aes(fact, log10(ratio)), fill=c("darkred","darkblue"), notch = T)+theme_bw()+ggtitle("lincRNA, 20%")+
  annotate(x=c(1, 2),y=c(-1,-1),geom = "text", label=c(short_med(whole_loc$ratio[whole_loc$fact=="loc_Tb_lincRNA20"]), 
                                                       short_med(whole_loc$ratio[whole_loc$fact=="loc_nTb_lincRNA20"])))+
  scale_x_discrete(name=paste("p-value = ", short_wilcox(whole_loc$ratio[whole_loc$fact=="loc_Tb_lincRNA20"],
                                                         whole_loc$ratio[whole_loc$fact=="loc_nTb_lincRNA20"]), sep=" "),
                   labels=c("TADb", "non-TADb"))
p5 <-ggplot(data=whole_loc[whole_loc$fact %in% c("loc_Tb_pc5","loc_nTb_pc5"),])+
  geom_boxplot(aes(fact, log10(ratio)), fill=c("darkred","darkblue"), notch = T)+theme_bw()+ggtitle("protein-coding, 5%")+
  annotate(x=c(1, 2),y=c(-1,-1),geom = "text", label=c(short_med(whole_loc$ratio[whole_loc$fact=="loc_Tb_pc5"]), 
                                                       short_med(whole_loc$ratio[whole_loc$fact=="loc_nTb_pc5"])))+
  scale_x_discrete(name=paste("p-value = ", short_wilcox(whole_loc$ratio[whole_loc$fact=="loc_Tb_pc5"],
                                                         whole_loc$ratio[whole_loc$fact=="loc_nTb_pc5"]), sep=" "),
                   labels=c("TADb", "non-TADb"))
p10 <-ggplot(data=whole_loc[whole_loc$fact %in% c("loc_Tb_pc10","loc_nTb_pc10"),])+
  geom_boxplot(aes(fact, log10(ratio)), fill=c("darkred","darkblue"), notch = T)+theme_bw()+ggtitle("protein-coding, 10%")+
  annotate(x=c(1, 2),y=c(-1,-1),geom = "text", label=c(short_med(whole_loc$ratio[whole_loc$fact=="loc_Tb_pc10"]), 
                                                       short_med(whole_loc$ratio[whole_loc$fact=="loc_nTb_pc10"])))+
  scale_x_discrete(name=paste("p-value = ", short_wilcox(whole_loc$ratio[whole_loc$fact=="loc_Tb_pc10"],
                                                         whole_loc$ratio[whole_loc$fact=="loc_nTb_pc10"]), sep=" "),
                   labels=c("TADb", "non-TADb"))
p20 <-ggplot(data=whole_loc[whole_loc$fact %in% c("loc_Tb_pc20","loc_nTb_pc20"),])+
  geom_boxplot(aes(fact, log10(ratio)), fill=c("darkred","darkblue"), notch = T)+theme_bw()+ggtitle("protein-coding, 20%")+
  annotate(x=c(1, 2),y=c(-1,-1),geom = "text", label=c(short_med(whole_loc$ratio[whole_loc$fact=="loc_Tb_pc20"]), 
                                                       short_med(whole_loc$ratio[whole_loc$fact=="loc_nTb_pc20"])))+
  scale_x_discrete(name=paste("p-value = ", short_wilcox(whole_loc$ratio[whole_loc$fact=="loc_Tb_pc20"],
                                                         whole_loc$ratio[whole_loc$fact=="loc_nTb_pc20"]), sep=" "),
                   labels=c("TADb", "non-TADb"))
grid.arrange(nrow=2, l5, l10, l20, p5, p10, p20)


