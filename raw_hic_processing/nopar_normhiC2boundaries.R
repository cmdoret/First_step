# In this script, I try to define TAD boundaries from normalized Hi-C data. I use the set of merged TADs to find borders and I set the boundaries
# limits where there is a decrease in interaction of...
# Cyril Matthey-Doret
# 20.10.2016

############################

# Loading data
#setwd("/Users/cmatthe5/Documents/First_step/data/")
setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
#setwd("/home/cyril/Documents/First_step/data/")

stopifnot(require(Matrix))
# Loading all matrices in a list (takes pretty long)
matlist <- list()
for(c in c("1")){
  tmp <- gzfile(paste0("norm_hic_data/chr",c,"_5kb_norm.txt.gz"),"rt")  # Files are unzipped before read
  matlist[[c]] <- c(readMM(tmp),c)
  close(tmp)  # Freeing connection for next matrix
}
TAD <- read.table("TAD/merged/merged_TAD.bed")
colnames(TAD) <- c("chr","start","end")

#==========================


#Defining functions

vec_sub_square <- function(v,s,e,n,w){  
  # this function allows to get the values inside a square sub matrix in a 1D vector representing a larger square matrix
  # v=vector,s=upper left corner of square, e=bottom right, n= number of cols in matrix,w=width of subsquare
  sub_sq <- c()
  c <- 1
  while(s<=e){
    for(i in s:(s+(w-1))){
      sub_sq <-append(sub_sq,v[i])
    }
    s <- s+n
  }
  return(sub_sq)
}

for_diam_slide<-function(m,R=5000, D=100000){  #Vectorized version of the slider. MUCH FASTER!
  L <- length(m[[1]][1,])
  M <- as.vector(t(m[[1]]))
  diam <- rep(0,L) # preallocating space for diamond-summed data.
  c <- 1
  for(r in (D/R+1):(L-D/R)){
    i = seq(from=(D/R+L*D/R+1),to=(L*L-(D/R+L*D/R)),by=(L+1))
    diam[r] <- sum(vec_sub_square(v=M,s=(i[c]-L*(D/R-1)),e=(i[c]+D/R-1),n=L,w=(D/R)))
    c <- c+1
  }
  return(diam)
}


vec_diam_slide<-function(m,R=5000, D=100000){  #Vectorized version of the slider. MUCH FASTER!
  L <- length(m[1,])
  M <- as.vector(t(m))
  diam <- rep(0,L) # preallocating space for diamond-summed data.
  diam[(D/R+1):(L-(D/R))] <- sapply(X = seq(from=(D/R+L*D/R+1),to=(L*L-(D/R+L*D/R)),by=(L+1)),
                                    simplify = T, FUN= function(d){
                                      sum(vec_sub_square(v=M,s=(d-L*(D/R-1)),e=(d+(D/R-1)),n=L,w=(D/R)))})
  return(diam)
}

find_bound<-function(m,R=5000, D=100000){ # D is diamond size, R is resolution
  sub_TAD <- TAD[as.character(TAD$chr)==paste0("chr",m[[2]]),]  # subsetting TAD dataframe by chromosome (1chr/node)
  sub_TAD$Lbound.start = sub_TAD$Lbound.end = sub_TAD$Rbound.start =
    sub_TAD$Rbound.end = sub_TAD$maxint <- rep(0,length(sub_TAD$start))  # allocating space for max interaction observed in each TAD.
  diam <- for_diam_slide(m[[1]],R,D)
  for(i in 1:length(sub_TAD$start)){  # Iterating over TADs
    start.ind <- pos_2_index(sub_TAD$start[i])  # transforming TAD start position to row in matrix. (TADs are already at 5kb res, so no need to approximate)
    end.ind <- pos_2_index(sub_TAD$end[i])  # Same for end position
    sub_TAD$maxint[i]<-max(diam[start.ind:end.ind])  # putting max interaction value for each TAD in dataframe.
    sl=el <- c(start.ind, diam[start.ind])  # Initiating left and right limits of left boundary
    while(sl[2]<(diam[start.ind]+sub_TAD$maxint[i]/10)){  # Start of left boundary -> sliding to the left
      sl[1] <- sl[1]-1; sl[2] <-diam[sl[1]]
      if(sl[1]==pos_2_index(D)){break}}
    sub_TAD$Lbound.start[i] <- index_2_pos(sl[1])  # + 1 because position corresponds to left edge of the bin. (don't want to include bin with interactions in boundary) ...forget it

    sub_TAD$Lbound.end[i] <- index_2_pos(el[1])
    sr=er <-c(end.ind, diam[end.ind])  # Initiating left and right limits of right boundary

    sub_TAD$Rbound.start[i] <- index_2_pos(sr[1])

    while(er[2]<(diam[end.ind]+sub_TAD$maxint[i]/10)){  # End of right boundary -> sliding to the right
      er[1] <- er[1]+1; er[2] <-diam[er[1]]
      if(er[1]==(length(m[[1]][1,])-(D/R))){break}}
    sub_TAD$Rbound.end[i] <- index_2_pos(er[1])
  }
  return(sub_TAD)
}

index_2_pos <- function(ind, resolution=5000) { return(resolution * (ind - 1)) }# Finding actual position in genome from bin number.
pos_2_index <- function(pos, resolution=5000) { return(1 + (pos / resolution))  }  # Finding index in vector from entry in matrix
#=========================


# Exporting the list of TADs and the required functions to each node in the cluster.

full <- find_bound(matlist[["1"]])
setwd("/home/cyril/Documents/Master/sem_1/First_step/raw_hic_processing/boundaries/")
write.table(full, file = "HICBOUND1.txt",quote = F,col.names = T,row.names = F,sep = "\t")

