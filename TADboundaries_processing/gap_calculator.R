# The purpose of this script is to compute all the gaps in a list of TADs and store them in a bed file. 
# This assumes TAD are not overlapping.
# Cyril Matthey-Doret
# 11.10.2016
############################

#Loading data
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



#This function generates a list of gaps from a bed file. The columns need to have specific names (chr, start, end)
gap_gen <- function(TADdf){
  gapdf <- data.frame()  # initiation of gap data frame
  for(c in levels(as.factor(TADdf$chr))){  # iterating over chromosomes
    chr_df <- TADdf[TADdf$chr==c,]   # Creating a subset of the dataframe for each chromosome
    chr_df <- chr_df[order(chr_df$start),] # Ordering the TADs according to their start site in each data frame
    for(r in row.names(chr_df)){  # Iterating over row names of TADs (assuming these are numbers and have not been changed)
      if(r < max(row.names(chr_df))){  # The last TAD of each chromosome has no gap after it.
        gapStart <- chr_df[r,"end"]  # Start of the gap is the end of the current TAD
        gapEnd <- chr_df[as.character(as.numeric(r)+1),"start"]  # End of the gap is the start of the next TAD
        gapdf <- rbind(gapdf,
                       cbind(paste("gap",r,as.character(as.numeric(r)+1),sep="_"),
                             c,
                             gapStart,
                             gapEnd))  
        # Each gap is added to the dataframe along with its chromosome number and an ID matching the surrounding TADs. For example, gap_23_24
        # is between TAD_23 and TAD_24
      }
    }
  }
  colnames(gapdf) <- c("ID","chr","start","end")
  return(gapdf)
}

# Calculating gaps
gaplist <- gap_gen(TAD)
write.table(gaplist,file = "TAD/merged/gaps.bed",sep="\t",quote = F,col.names = F,row.names = F)
