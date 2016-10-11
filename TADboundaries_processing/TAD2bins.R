# The purpose of this script is to divide TADs into bins 
#Cyril Matthey-Doret
#11.10.2016

setwd("/home/cyril/Documents/Master/sem_1/First_step/data")
#Loading domains and boundaries.
TADb10 <- read.table("TAD/merged/TAD_boundaries10.bed")
#TADb10 <- TADb10[,-4]
TAD <- read.table("TAD/merged/merged_TAD.bed")
TAD <- cbind(paste("merged",row.names(TAD),sep="_"),TAD)
#TAD <- TAD[,-4]
colnames(TADb10) <-c("chr", "start","end")
colnames(TAD)<-c("ID","chr", "start","end")
# Loading lincRNAs and protein coding genes.
nTADb_lincRNA10 <- read.table("linc_RNA/merged/nonTADbound-lincRNA10.bed")
TADb_lincRNA10 <- read.table("linc_RNA/merged/nonTADbound-lincRNA10.bed")
TADb_pcgene10 <- read.table("pc_genes/TADbound-pcgene10.bed")
nTADb_pcgene10 <- read.table("pc_genes/nonTADbound-pcgene10.bed")
colnames(nTADb_lincRNA10)=colnames(TADb_pcgene10)=
  colnames(TADb_lincRNA10) = colnames(nTADb_pcgene10)<-c("chr", "start","end","ID", "strand")
gaps <- read.table("TAD/merged/gaps.bed")
colnames(gaps) <- c("ID", "chr", "start", "end")

#==========================================

# Division of TADs: each TAD must be divided into 10 bins, plus 3 bins outside.
# The length of each bin should normally be 10% of the TAD length, but for outer bins, 
# We will use the smallest value between 3*10% and 0.5*gap
TAD_bins <- data.frame()
TAD.splitter<-function(tad,gap){ #the function takes a list of TADs and the corresponding list of gaps.
  tad_nr<-strsplit(as.character(tad[1]),split = "_",fixed = T)[[1]][2]
  # Splitting id of the gap (e.g. gap_3_4). gap with first number = tad will be after, gap with second number = tad will be before
  next.gap <- gap[as.numeric(strsplit(as.character(gap[,1]),split = "_")[[1]][2])==as.numeric(tad_nr),] 
  prev.gap <- gap[as.numeric(strsplit(as.character(gap[,1]),split = "_")[[1]][3])==as.numeric(tad_nr),]
  # Initiating left outer bins.
  if((prev.gap/2)<(tad$end-tad$start)*0.05){
    #no outer bins
  }
  else if((prev.gap/2)<(tad$end-tad$start)*0.1){
    #only 0.05 outer bins
  }else if((prev.gap/2)<(tad$end-tad$start)*0.2){
    #only 0.1 and 0.05 outer bins
  }
  else{
    #0.2,0.1 and 0.05 outer bins
  }
  # Initiating right outer bins
  if((next.gap/2)<(tad$end-tad$start)*0.05){
    #no outer bins
  }
  else if((next.gap/2)<(tad$end-tad$start)*0.1){
    #only 0.05 outer bins
  }else if((next.gap/2)<(tad$end-tad$start)*0.2){
    #only 0.1 and 0.05 outer bins
  }
  else{
    #0.2,0.1 and 0.05 outer bins
  }
  #initiating inner bins
  for(i in seq(from=0,to=1,by=0.1)){
    # split into 10 bins
  }
}

apply(TAD, MARGIN = 1,FUN=TAD.splitter, gap=gaps)
