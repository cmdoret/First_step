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

find_bound<-function(m){
  sub_TAD <- TAD[as.character(TAD$chr)==paste0("chr",m[2]),]
  sub_TAD$maxint <- rep(0,length(sub_TAD$start))
  for(i in 1:length(sub_TAD$start)){
    start.ind <- pos_2_index(sub_TAD$start[i])
    end.ind <- pos_2_index(sub_TAD$end[i])
    sub_TAD$maxint[i]<-max(rowSums(m[1][[1]][start.ind:end.ind,]))
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
clusterExport(clus, varlist=c("TAD", "index_2_pos", "pos_2_index"), envir=environment())  # Exporting the list of TADs to each node in the cluster.
tmp <- clusterApplyLB( clus, matlist, find_bound)
