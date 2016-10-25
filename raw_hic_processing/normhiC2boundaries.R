# In this script, I try to define TAD boundaries from normalized Hi-C data. I use the set of merged TADs to find borders and I set the boundaries
# limits where there is a decrease in interaction of...
# Cyril Matthey-Doret
# 20.10.2016

############################

# Loading data
#setwd("/Users/cmatthe5/Documents/First_step/data/")
setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
#setwd("/home/cyril/Documents/First_step/data/")
if (!exists("n.cores")) { 
  stopifnot(require(parallel)) 
  n.cores <- parallel:::detectCores()-2 # Value depends on number of cores in the CPU
}
stopifnot(require(snow))
stopifnot(require(Matrix))
# Loading all matrices in a list (takes pretty long)
matlist <- list()
for(c in c(22)){
  tmp <- gzfile(paste0("norm_hic_data/chr",c,"_5kb_norm.txt.gz"),"rt")  # Files are unzipped before read
  matlist[[c]] <- c(readMM(tmp),c)
  close(tmp)  # Freeing connection for next matrix
}
TAD <- read.table("TAD/merged/merged_TAD.bed")
colnames(TAD) <- c("chr","start","end")

#==========================

clus <- makeSOCKcluster( n.cores  )  # Creating cluster
on.exit( stopCluster(clus) )  

#Defining functions

find_bound<-function(m,resolution=5000, diam_size=100000){
  sub_TAD <- TAD[as.character(TAD$chr)==paste0("chr",m[[2]]),]  # subsetting TAD dataframe by chromosome (1chr/node)
  sub_TAD$Lbound.start = sub_TAD$Lbound.end = sub_TAD$Rbound.start =
    sub_TAD$Rbound.end = sub_TAD$maxint <- rep(0,length(sub_TAD$start))  # allocating space for max interaction observed in each TAD.
  diam <- data.frame(ind=seq(1,length(m[[1]][1,])),int=rep(0,length(m[[1]][1,])))  # preallocating space for diamond-summed data.
  #for(d in 1:length(m[[1]][1,])){diam$int[d] <- sum(rowSums(m[[1]][d:(d+(diam_size/resolution)),d:(d+(diam_size/resolution))]))}
  diam$int[pos_2_index(diam_size):(length(m[[1]][1,])-(diam_size/resolution))] <- sapply(X = seq(pos_2_index(diam_size),
                                                                                                length(m[[1]][1,]))-(diam_size/resolution),
                                                                                        simplify = T, FUN= function(d){
    sum(rowSums(m[[1]][d:(d+(diam_size/resolution)),d:(d+(diam_size/resolution))]))})
  for(i in 1:length(sub_TAD$start)){  # Iterating over TADs
    start.ind <- pos_2_index(sub_TAD$start[i])  # transforming TAD start position to row in matrix. (TADs are already at 5kb res, so no need to approximate)
    end.ind <- pos_2_index(sub_TAD$end[i])  # Same for end position
    sub_TAD$maxint[i]<-max(diam$int[start.ind:end.ind])  # putting max interaction value for each TAD in dataframe.
    sl=el <- c(start.ind, diam$int[start.ind])  # Initiating left and right limits of left boundary
    while(sl[2]<sub_TAD$maxint[i]/10){  # Start of left boundary -> sliding to the left
      sl[1] <- sl[1]-1; sl[2] <-diam$int[sl[1]]}
    sub_TAD$Lbound.start[i] <- index_2_pos(sl[1])  # + 1 because position corresponds to left edge of the bin. (don't want to include bin with interactions in boundary) ...forget it

    while(el[2]<sub_TAD$maxint[i]/10){  # End of left boundary -> sliding to the right
      el[1] <- el[1]+1; el[2] <-diam$int[el[1]]}
    sub_TAD$Lbound.end[i] <- index_2_pos(el[1])
    sr=er <-c(end.ind, diam$int[end.ind])  # Initiating left and right limits of right boundary

    while(sr[2]<sub_TAD$maxint[i]/10){  # Start of right boundary -> sliding to the left
      sr[1] <- sr[1]-1; sr[2] <-diam$int[sr[1]]}
    sub_TAD$Rbound.start[i] <- index_2_pos(sr[1])

    while(er[2]<sub_TAD$maxint[i]/10){  # End of right boundary -> sliding to the right
      er[1] <- er[1]+1; er[2] <-diam$int[er[1]]}
    sub_TAD$Rbound.end[i] <- index_2_pos(er[1])
  }
  return(sub_TAD)
}

index_2_pos <- function(ind, resolution=5000) { return(resolution * (ind - 1)) }# Finding actual position in genome from bin number.
pos_2_index <- function(pos, resolution=5000) { return(1 + (pos / resolution))  }  # Finding index in vector from entry in matrix
#=========================

# Sending tasks
clusterCall( clus, function() { 
  require(Matrix)       # executed on each node. 
  # Loading required libraries in each instance of R
} )
# Exporting the list of TADs and the required functions to each node in the cluster.
clusterExport(clus, c("TAD", "index_2_pos", "pos_2_index"), envir=environment())
tmp <- clusterApplyLB( clus, matlist, find_bound)
full <- do.call("rbind", tmp)
full
write(full, file="BOUNDARIES_HIC.txt", sep= "\t")
#write.table(full, file = "TAD/merged/hic_bound.txt",quote = F,col.names = F,row.names = F,sep = "\t")