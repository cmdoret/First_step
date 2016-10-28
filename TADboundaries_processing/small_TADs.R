# This script processes a set of TADs to remove the large encompassing ones and keep only those that do not completely overlap.
# Cyril Matthey-Doret
# Fri Oct 28 08:49:55 2016 ------------------------------

# Loading data:
setwd("/home/cyril/Documents/First_step/data/")
setwd("/Users/cmatthe5/Documents/First_step/data/")
TAD <- read.table("TAD/GM12878_TAD_domains.bed")
TAD <- TAD[,c(1,2,3,5)]
colnames(TAD) <- c("chr","start","end","ID")
library(intervals)
reduce_TAD <- function(m){
  blacklist <- c()
  chr<-TAD$chr
  start <- TAD$start
  end<-TAD$end
  id <- as.character(TAD$ID)
  for (b in 1:length(TAD$ID)){
    if(all(id[b]!=m[4], 
           chr[b]==m[1],
           start[b]<=m[2],
           end[b]>=m[3])){
      blacklist <- append(blacklist,as.character(TAD$ID[b]))
    }
  }
  return(blacklist)
}
long_TAD <- apply(TAD,MARGIN = 1,FUN = reduce_TAD)
#======================================
# with overlap package

reduce_TAD <- function(TAD){
  blacklist <- list()
  for (c in TAD$chr){  # Doing operations separately by chromosomes
    blacklist[[c]] <- c()
    TAD_tmp<-Intervals(matrix(c(TAD$start[TAD$chr==c],
                       TAD$end[TAD$chr==c]),ncol=2))
    rownames(TAD_tmp) <- TAD$ID[TAD$chr==c]  # Naming intervals after TADs
    over <-interval_overlap(TAD_tmp,TAD_tmp)  # Calling overlaps between all TADs
    for(k in names(over)){  # Iterating over TAD names
      if(length(over[[k]])>1){  # If the number of overlaps for a TAD is gt 1 (i.e. overlaps not only with itself)
        for(t in over[[k]]){  # Iterating over items overlapping this TAD
          if(k!=names(over[t]) & size(TAD_tmp[k]<=size(TAD_tmp[t]))){  # Only counting overlaps that are not the TAD itself
            blacklist[[c]] <- append(blacklist[[c]],rownames(TAD_tmp[t]))
          }
        }
      }
    }
  }
  return(blacklist)
}

large <-reduce_TAD(TAD)



