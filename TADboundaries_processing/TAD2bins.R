# The purpose of this script is to divide TADs into bins 
#Cyril Matthey-Doret
#11.10.2016

#setwd("/home/cyril/Documents/First_step/data/")
setwd("/home/cyril/Documents/Master/sem_1/First_step/data")
#setwd("/Users/cmatthe5/Documents/First_step/data/")
#Loading domains and boundaries.
TADb10 <- read.table("TAD/merged/TAD_boundaries10.bed")
#TADb10 <- TADb10[,-4]
TAD <- read.table("TAD/short/short_TADs.bed")
TAD <- cbind(paste("short",row.names(TAD),sep="_"),TAD)
TAD <- TAD[,-5]
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
# We will use the smallest value between 3*10% and 3*gap/6
TAD_bins <- data.frame()
TAD.splitter<-function(tad,gap,bw){ #the function takes a list of TADs and the corresponding list of gaps.
  start_tad <- as.numeric(unname(tad[3]))
  end_tad<-as.numeric(unname(tad[4]))
  tad_nr<-strsplit(as.character(tad[1]),split = "_",fixed = T)[[1]][2]
  # Splitting id of the gap (e.g. gap_3_4). gap with first number = tad will be after, gap with second number = tad will be before
  gap_nrs<-strsplit(as.character(gap[,1]),split = "_")
  idx<-1
  prev.gap <- NA
  next.gap <- NA
  for(g in gap_nrs){
      if(as.numeric(g[3])==as.numeric(tad_nr)){prev.gap <- gap[idx,]}
      if(as.numeric(g[2])==as.numeric(tad_nr)){next.gap <- gap[idx,]}
    #print(head(gap_nrs))
    idx <- idx+1
  }
  # bins will take the form: c(TAD_ID,chr,bin_nr, start, end)
  thr <- seq(from=bw,to=0.3,by=bw)
  Lbin<-data.frame()
  Rbin<-data.frame()
  for(t in 1:length(thr)){
    Lbin <- rbind(Lbin, rep(NA,5))
    Rbin <- rbind(Rbin, rep(NA,5))
  }
  prev.gap_len <-as.numeric(prev.gap[4])-as.numeric(prev.gap[3])
  next.gap_len <-as.numeric(next.gap[4])-as.numeric(next.gap[3]) 
  # Initiating left outer bins.
  if(all(!is.na(prev.gap), tad[2]==paste0("chr",prev.gap[2]))){  # Checking if TAD and left gap are on same chromosome
    if((prev.gap_len/2)>(end_tad-start_tad)*0.3){
      c <- 1
      for(t in thr){
        Lbin[c,] <-c(tad[1],tad[2],paste0("L",c),start_tad-(end_tad-start_tad)*t,start_tad-(end_tad-start_tad)*(t-thr[1]))
        c <- c+1
      }
    }
    else{  # Half the gap is smaller than 30% of TAD length -> define outer bins as 1/6th of the gap each.
      c <- 2
      gap_thr <- ceiling(seq(from=0,to=prev.gap_len/2,by=(prev.gap_len/2)/length(thr)))
      for(t in gap_thr[2:length(gap_thr)]){
        Lbin[c-1,] <-c(tad[1],tad[2],paste0("L",(c-1)),start_tad-t,start_tad-gap_thr[c-1])
        c <- c+1
      }
    }
  }
  else{  # TAD and gap on separate chromosomes; all 3 left outer bins can be used
    c <- 1
    for(t in thr){
      Lbin[c,] <-c(tad[1],tad[2],paste0("L",c),start_tad-(end_tad-start_tad)*t,start_tad-(end_tad-start_tad)*(t-thr[1]))
      c <- c+1
    }
  }
  # Initiating right outer bins
  if(all(!is.na(next.gap), tad[2]==paste0("chr",next.gap[2]))){   # Checking if TAD and right gap are on same chromosome
    if((next.gap_len/2)>(end_tad-start_tad)*0.3){
      c <- 1
      for(t in thr){
        Rbin[c,] <-c(tad[1],tad[2],paste0("R",c),end_tad+(end_tad-start_tad)*(t-thr[1]),end_tad+(end_tad-start_tad)*t)
        c <- c+1
      } 
    }
    else{  # Half the gap is smaller than 30% of TAD length -> define outer bins as 1/6th of the gap each.
      c <- 2
      gap_thr <- floor(seq(from=0,to=next.gap_len/2,by=(next.gap_len/2)/length(thr)))
      for(t in gap_thr[2:length(gap_thr)]){
        Rbin[c-1,] <-c(tad[1],tad[2],paste0("R",(c-1)),end_tad+gap_thr[c-1],end_tad+t)
        c <- c+1
      }
    }
  }
  else{   #TAD and gap on separate chromosomes; all 3 left outer bins can be used
    c <- 1
    for(t in thr){
      Rbin[c,] <-c(tad[1],tad[2],paste0("R",c),end_tad+(end_tad-start_tad)*(t-thr[1]),end_tad+(end_tad-start_tad)*t)
      c <- c+1
    }
  }
  ibin <- data.frame()  #initiating inner bins
  for(b in 1:1/thr[1]){
    ibin <-rbind(ibin,c(rep(NA,5)))
  }
  c <- 1
  for(i in seq(from=thr[1],to=1,by=thr[1])){
    # split into 10 bins.
    ibin[c,] <- c(tad[1],tad[2],as.character(i*100),start_tad+(end_tad-start_tad)*(i-thr[1]),start_tad+(end_tad-start_tad)*(i))
    c <- c+1
  }
  # appending all bins for a given TAD to the dataframe
  return(list(Lbin,Rbin,ibin))
}

TAD_bins <-apply(TAD, MARGIN = 1,FUN=TAD.splitter, gap=gaps,bw=0.01)
#bw parameter sets the binwidth
#=========================================================================
#Only inner bins:

inner_bins <- function(TAD,w){  # Takes a TAD as an input and the binwidth as %TAD length (e.g. 0.05 for 5%)
  TAD_start <- as.numeric(TAD[3])
  TAD_end <-   as.numeric(TAD[4])
  TAD_ID <- as.character(TAD[1])
  TAD_chr <- as.character(TAD[2])
  bin_coord <- seq(TAD_start,TAD_end,(TAD_end - TAD_start)*w)
  nbins <- length(bin_coord)-1
  bins <- data.frame(ID=rep(NA,nbins),
                     chr=rep(NA,nbins),
                     start=rep(NA,nbins),
                     end=rep(NA,nbins),
                     bin=rep(NA,nbins))
  for(i in 1:nbins){
    bins[i,] <- c(TAD_ID,TAD_chr,bin_coord[i],bin_coord[i+1],i)
  }
  return(bins)  # Returns all inner bins for the input TAD
}
options(scipen=999)
short_bins <-apply(X=TAD,MARGIN = 1,FUN=inner_bins, w=0.05)
short_bins <- do.call("rbind",short_bins)
write.table(short_bins,file="TAD/short/short_bins5.bed",quote = F, sep="\t",col.names = F,row.names = F)
#========================================================================

#Concatenation of all lists into a single humongus dataframe
k = length(unlist(TAD_bins))/5
whole_bins <- data.frame(rep(0,times=k),rep(0,times=k),rep(0,times=k),rep(0,times=k),rep(0,times=k))
library(data.table)
c <-1
v <-length(unlist(TAD_bins[[1]]))/5
for(i in TAD_bins){
  whole_bins[c:v,] <- as.data.frame(rbindlist(i))
  c <- c+length(unlist(TAD_bins[[1]]))/5
  v <- v+length(unlist(TAD_bins[[1]]))/5
}

options(scipen=999)
write.table(na.omit(whole_bins), quote=F,sep="\t",row.names = F,col.names = F,file = "TAD/merged/merged_TADbins1.txt")
options(scipen=0)
#===============================================================================
# formatting file for bedtools
whole_bins <- read.table("TAD/merged/merged_TADbins10.txt")
whole_bins <- whole_bins[,c(2,4,5,1,3)]
write.table(whole_bins, quote=F,sep="\t",row.names = F,col.names = F,file = "TAD/merged/short_TADbins10.bed")
