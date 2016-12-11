# This script processes a set of TADs to remove the large encompassing ones and keep only those that do not completely overlap.
# Cyril Matthey-Doret
# Fri Oct 28 08:49:55 2016 ------------------------------

# Loading data:
setwd("/home/cyril/Documents/First_step/data/")
setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
setwd("/Users/cmatthe5/Documents/First_step/data/")
TAD <- read.table("TAD/GM12878_TAD_domains.bed")
TAD <- TAD[,c(1,2,3,5)]
colnames(TAD) <- c("chr","start","end","ID")
library(intervals)

##############################

# Only discarding TADs if a smaller one is entirely inside.
inbig_TAD <- function(TAD){
  blacklist <- list()
  for (c in levels(TAD$chr)){  # Doing operations separately by chromosomes
    blacklist[[c]] <- c()
    TAD_tmp<-Intervals(matrix(c(TAD$start[TAD$chr==c],
                                TAD$end[TAD$chr==c]),ncol=2))
    rownames(TAD_tmp) <- TAD$ID[TAD$chr==c]  # Naming intervals after TADs
    over <-interval_overlap(TAD_tmp,TAD_tmp)  # Calling overlaps between all TADs
    for(k in names(over)){  # Iterating over TAD names
      if(length(over[[k]])>1){  # If the number of overlaps for a TAD is gt 1 (i.e. overlaps not only with itself)
        for(t in over[[k]]){  # Iterating over items overlapping this TAD
          if(k!=names(over[t]) # Only counting overlaps that are not the TAD itself
             & size(TAD_tmp[k])<=size(TAD_tmp[t])   # Removing TAD if larger or same size.
             & size(interval_intersection(TAD_tmp[k],TAD_tmp[t]))>=size(TAD_tmp[k])){  # Removing TAD if it completely encompass smaller one
            blacklist[[c]] <- append(blacklist[[c]],rownames(TAD_tmp[t])) # Contains only large TADs that need to be removed
          }
        }
      }
    }
  }
  return(blacklist)
}
large <-inbig_TAD(TAD)  # Storing the IDs of all large TADs

options(scipen=999)  # Preventing R from using scientific notation (other programs such as bedtools and genome browsers do not like it)
short <- TAD[!(TAD$ID %in% unname(unlist(large))),]  # Removing all large TADs from the original list
write.table(short,"TAD/short_TADs.bed",sep="\t",quote = F,col.names=F,row.names=F)


