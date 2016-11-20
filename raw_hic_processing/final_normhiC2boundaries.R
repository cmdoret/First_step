# In this script, I try to define TAD boundaries from normalized Hi-C data. I use the set of merged TADs to find borders 
# and I set the boundaries limits based on changes in interactions.
# This version of the scripts splits the calculation of interactions and the boundaries definition into 2 separate tasks to make it
# less prone to crashes.
# Cyril Matthey-Doret
# 20.10.2016

############################

# Loading data
#setwd("/Users/cmatthe5/Documents/First_step/data/")
#setwd("/scratch/beegfs/monthly/mls_2016/cmatthey/first_step/data/")
#setwd("/home/cyril/Documents/First_step/data/")
setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")

stopifnot(require(Matrix))
# Loading all matrices in a list (takes pretty long)
matlist <- list()
for(c in c("3")){
  tmp <- gzfile(paste0("norm_hic_data/GM12878/chr",c,"_5kb_norm.txt.gz"),"rt")  # Files are unzipped before read
  matlist[[c]] <- c(readMM(tmp),c)
  close(tmp)  # Freeing connection for next matrix
}
TAD <- read.table("TAD/short_fullover_TAD.bed")
colnames(TAD) <- c("chr","start","end")

#==========================


#Defining functions

D2toD1 <- function(I,J,N){
  return(((I-1)*N)+J)
}

vec_sub_square <- function(v,s,e,n,w){  
  # this function allows to get the values inside a square sub matrix in a 1D vector representing a larger square matrix
  # v=vector,s=upper left corner of square, e=bottom right, n= number of cols in matrix,w=width of subsquare
  sub_sq <- rep(0,times=w*w)
  cs <- 1
  while(s<=e){
    for(i in s:(s+(w-1))){
      sub_sq[cs] <-v[i]
      cs <- cs+1
    }
    s <- s+n
  }
  return(sub_sq)
}

vec_diam_slide<-function(m,R=5000, D=100000){  #Vectorized version of the slider. MUCH FASTER!
  L <- length(m[1,])
  M <- as.vector(t(m))
  diam <- rep(0,L) # preallocating space for diamond-summed data.
  diam[(D/R):(L-(D/R-1))] <- sapply(X = seq(from=D2toD1(D/R,D/R,L),to=D2toD1((L-(D/R-1)),(L-(D/R-1)),L),by=(L+1)),
                                    simplify = T, FUN= function(d){
                                      sum(vec_sub_square(v=M,s=(d-L*(D/R-1)),e=(d+(D/R-1)),n=L,w=(D/R)))})
  return(diam)
}

find_bound<-function(n,R=5000, D=100000){ # D is diamond size, R is resolution
  sub_TAD <- TAD[as.character(TAD$chr)==paste0("chr",n),]  # subsetting TAD dataframe by chromosome (1chr/node)
  print(n)
  sub_TAD$Lbound.start = sub_TAD$Lbound.end = sub_TAD$Rbound.start =
    sub_TAD$Rbound.end = sub_TAD$maxint <- rep(0,length(sub_TAD$start))  # allocating space for max interaction observed in each TAD.
  diam <- scan(file = paste0("diam_sums/GM12878/chr",n),what = numeric())
  for(i in 1:length(sub_TAD$start)){  # Iterating over TADs
    start.ind <- pos_2_index(sub_TAD$start[i])  # transforming TAD start position to row in matrix. (TADs are already at 5kb res, so no need to approximate)
    end.ind <- pos_2_index(sub_TAD$end[i])  # Same for end position
    sub_TAD$maxint[i]<-max(diam[start.ind:end.ind])  # putting max interaction value for each TAD in dataframe.
    sl=el <- c(start.ind, diam[start.ind])  # Initiating left and right limits of left boundary

    sub_TAD$Lbound.start[i] <- index_2_pos(sl[1])  # + 1 because position corresponds to left edge of the bin. (don't want to include bin with interactions in boundary) ...forget it
    while(el[2]<(diam[start.ind]+sub_TAD$maxint[i]/10)){  # End of left boundary -> sliding to the right
      el[1] <- el[1]+1; el[2] <-diam[el[1]]
      if(el[1]==(length(diam)-(D/R))){break}}
    sub_TAD$Lbound.end[i] <- index_2_pos(el[1])
    sr=er <-c(end.ind, diam[end.ind])  # Initiating left and right limits of right boundary
    
    while(sr[2]<(diam[end.ind]+sub_TAD$maxint[i]/10)){  # Start of right boundary -> sliding to the left
      sr[1] <- sr[1]-1; sr[2] <-diam[sr[1]]
      if(sr[1]==pos_2_index(D)){break}}
    sub_TAD$Rbound.start[i] <- index_2_pos(sr[1])

    sub_TAD$Rbound.end[i] <- index_2_pos(er[1])
  }
  return(sub_TAD)
}

index_2_pos <- function(ind, resolution=5000) { return(resolution * (ind - 1)) }# Finding actual position in genome from bin number.
pos_2_index <- function(pos, resolution=5000) { return(1 + (pos / resolution))  }  # Finding index in vector from entry in matrix
#=========================

# 2 steps method:

# STEP 1: Only calculating the sums of interactions in each window and storing them: 
# If it still crashes, see alternative STEP 1 at the end.
for(n in c("9")){
  i <- matlist[[n]][[1]]
  diam<- vec_diam_slide(i)
  write.table(diam, file=paste0("diam_sums/GM12878/chr",matlist[[n]][[2]]),sep = "\t",quote=F,row.names = F,col.names = F)
}

# STEP 2: loading summed interactions and defining boundaries from the vectors.
for(n in c(seq(21,22))){
  bound<- find_bound(n)
  #write.table(diam, file=paste0("diam_sums/GM12878/chr",matlist[[n]][[2]]),sep = "\t",quote=F,row.names = F,col.names = F)
  print(bound)
}


##################################################
# Alternative STEP 1 in case the matrix is too large and needs to be splitted into overlapping parts:
# NOT TESTED!

D <- 100000; R <- 5000; P <- D/R # The size of the sliding diamond and the resolution will determine how the 2 matrices should overlap.

for(n in c("3")){
  k <- dim(matlist[["3"]][[1]])[1]
  i <- matlist[[n]][[1]][1:(k/2 + P),1:(k/2 + P)]
  diam.i<- vec_diam_slide(i)
  j <- matlist[[n]][[1]][(k/2 - P):k,(k/2 - P):k]
  diam.j<- vec_diam_slide(j)
  diam <- append(diam.i[-((length(diam.i)-(P-1)):length(diam.i))],diam.j[-(1:(P+1))])
  write.table(diam, file=paste0("diam_sums/GM12878/chr",matlist[[n]][[2]]),sep = "\t",quote=F,row.names = F,col.names = F)
}
