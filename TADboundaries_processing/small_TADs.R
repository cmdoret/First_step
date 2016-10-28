# This script processes a set of TADs to remove the large encompassing ones and keep only those that do not completely overlap.
# Cyril Matthey-Doret
# Fri Oct 28 08:49:55 2016 ------------------------------

# Loading data:
setwd("/home/cyril/Documents/First_step/data/")
TAD <- read.table("TAD/GM12878_TAD_domains.bed")
TAD <- TAD[,c(1,2,3,5)]
colnames(TAD) <- c("chr","start","end","ID")
blacklist <- c()
for(t in 1:length(TAD$ID)){
  for (b in 1:length(TAD$ID)){
    if(all(TAD[b,"ID"]!=TAD[t,"ID"], 
           TAD[b,"start"]<=TAD[t,"start"],
           TAD[b,"end"]>=TAD[t,"end"],
           TAD[b,"chr"]==TAD[t,"chr"])){
      blacklist <- append(blacklist,TAD$ID[b])
      }
  }
}
short_TAD <- TAD[!(TAD$ID %in% blacklist)]

blacklist <- c()
for(t in 1:length(TAD$ID)){
  for (b in 1:length(TAD$ID)){
    if(TAD[b,"ID"]!=TAD[t,"ID"] &
           TAD[b,"start"]<=TAD[t,"start"] &
           TAD[b,"end"]>=TAD[t,"end"] &
           TAD[b,"chr"]==TAD[t,"chr"]){
      blacklist <- append(blacklist,TAD$ID[b])
    }
  }
}
short_TAD <- TAD[!(TAD$ID %in% blacklist)]