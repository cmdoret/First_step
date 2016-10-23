# In this script, I try to define TAD boundaries from normalized Hi-C data. I use the set of merged TADs to find borders and I set the boundaries
# limits where there is a decrease in interaction of...
# Cyril Matthey-Doret
# 20.10.2016
#BSUB -o hicbound_out.txt
#BSUB -e hicbound_errors.txt
#BSUB -J hic_boundaries

############################

# Loading data
#setwd("/Users/cmatthe5/Documents/First_step/data/")
setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")

if (!exists("n.cores")) { 
  stopifnot(require(parallel)) 
  n.cores <- parallel:::detectCores() # Value depends on number of cores in the CPU
}
stopifnot(require(snow))
stopifnot(require(Matrix))
# Loading all matrices in a list (takes pretty long)
matlist <- list()
for(c in c(22,"X")){
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

find_bound<-function(m){
  sub_TAD <- TAD[as.character(TAD$chr)==paste0("chr",m[2]),]  # subsetting TAD dataframe by chromosome (1chr/node)
  sub_TAD$maxint <- rep(0,length(sub_TAD$start))  # allocating space for max interaction observed in each TAD.
  for(i in 1:length(sub_TAD$start)){  # Iterating over TADs
    start.ind <- pos_2_index(sub_TAD$start[i])  # transforming TAD start position to row in matrix. (TADs are already at 5kb res, so no need to approximate)
    end.ind <- pos_2_index(sub_TAD$end[i])  # Same for end position
    sub_TAD$maxint[i]<-max(rowSums(m[1][[1]][start.ind:end.ind,]))  # putting max interaction value for each TAD in dataframe.
    # compute threshold based on interactions at border and max interactions in TAD
    # slide before and after border to find bins corresponding to threshold
    # corresponding bins will be the limits of TAD boundaries. Keep in mind the coordinate obtained with index_2_pos(i) is bin's left edge.
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
