# The purpose of this script is to divide TADs into bins 
#Cyril Matthey-Doret
#11.10.2016

#setwd("/home/cyril/Documents/First_step/data/")
setwd("/home/cyril/Documents/Master/sem_1/First_step/data")
#setwd("/Users/cmatthe5/Documents/First_step/data/")
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
