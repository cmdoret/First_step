# The purpose of this script is to divide TADs into equally sized bins. The default length of bins is 10% of TAD length but 
# it can be changed by modifying the "w" parameter.
# Cyril Matthey-Doret
# 11.10.2016

setwd("/home/cyril/Documents/Master/sem_1/First_step/data")
#Loading domains 
TAD <- read.table("TAD/short/short_fullover_TAD.bed")
TAD <- cbind(paste("short",row.names(TAD),sep="_"),TAD)
TAD <- TAD[,-5]
colnames(TAD)<-c("ID","chr", "start","end")
# Loading lincRNAs and protein coding genes.


#==========================================
#Only inner bins:

inner_bins <- function(TAD,w){  # Takes a TAD as an input and the binwidth as %TAD length (e.g. 0.05 for 5%)
  TAD_start <- as.numeric(TAD[3])  # Making sure all values are in the right type.
  TAD_end <-   as.numeric(TAD[4])
  TAD_ID <- as.character(TAD[1])
  TAD_chr <- as.character(TAD[2])
  bin_coord <- seq(TAD_start,TAD_end,(TAD_end - TAD_start)*w)  #  Storing the start/end coordinates of each bin
  nbins <- length(bin_coord)-1  # Getting the number of bins
  bins <- data.frame(ID=rep(NA,nbins),  # Preallocating a dataframe for storing bins. There are as many rows as there are bins.
                     chr=rep(NA,nbins),
                     start=rep(NA,nbins),
                     end=rep(NA,nbins),
                     bin=rep(NA,nbins))
  for(i in 1:nbins){  # iterating on rows of the new dataframe
    bins[i,] <- c(TAD_ID,TAD_chr,bin_coord[i],bin_coord[i+1],i)  # Filling each row with the coordinates of the corresponding bin 
  }
  return(bins)  # Returns all inner bins for the input TAD
}
options(scipen=999)
short_bins <-apply(X=TAD,MARGIN = 1,FUN=inner_bins, w=0.1)  # Applying the innerbins function to all TADs in the list
short_bins <- do.call("rbind",short_bins)  # Concatenating results into a single dataframe.
write.table(short_bins,file="TAD/short/short_fullover_bins10.bed",quote = F, sep="\t",col.names = F,row.names = F)

