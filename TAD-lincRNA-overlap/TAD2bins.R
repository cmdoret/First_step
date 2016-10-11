# The purpose of this script is to divide TADs into bins 
setwd("/home/cyril/Documents/First_step/data")
#Loading domains and boundaries.
TADb10 <- read.table("TAD/merged/TAD_boundaries10.bed")
#TADb10 <- TADb10[,-4]
TAD <- read.table("TAD/merged/merged_TAD.bed")
#TAD <- TAD[,-4]
colnames(TADb10) = colnames(TAD)<-c("chr", "start","end")
# Loading lincRNAs and protein coding genes.
nTADb_lincRNA10 <- read.table("linc_RNA/merged/nonTADbound-lincRNA10.bed")
TADb_lincRNA10 <- read.table("linc_RNA/merged/nonTADbound-lincRNA10.bed")
TADb_pcgene10 <- read.table("pc_genes/TADbound-pcgene10.bed")
nTADb_pcgene10 <- read.table("pc_genes/nonTADbound-pcgene10.bed")
colnames(nTADb_lincRNA10)=colnames(TADb_pcgene10)=
  colnames(TADb_lincRNA10) = colnames(nTADb_pcgene10)<-c("chr", "start","end","ID", "strand")

#======================================================================
# Calculating gaps

TADdf <- TAD[order(TAD$chr,TAD$start),]
gaplist <- by(TADdf,"chr",gap_gen)


gap_gen <- function(TADdf){
  gaps<-c()
  for(i in row.names(TADdf)){
    gapstart <- TADdf$end
    gaps <- append(gaps,gap)
  }
  return(gaps)
}
