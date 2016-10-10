# This script intends to cut TAD boundaries so that a boundary does not overlap with the boundary of the next TAD.

# Cyril Matthey-Doret 
# 07.10.2016
######################################

setwd("/home/cyril/Documents/First_step/data/")

TADb10 <- read.table("TAD/TAD_boundaries10.bed")
TADb10 <- TADb10[,-4]
colnames(TADb10)<-c("chr", "start","end","ID")

#by(chr,apply(TADb,fun))

over.bound <- function(b1,b2){
  if(any(b1$start > b2$start & b1$end < b2$end, 
         b1$start < b2$start & b1$end > b2$end)){   # If one boundary is entirely contained in the other.
    
  }
}
  
cut.bound<-function(df){
  for(id in df$ID){
    for(o in df$ID[-(df$ID==id),]){
      
    }
  }
}