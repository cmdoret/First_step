# In this script, I define TAD boundaries from normalized Hi-C data. I use the set of filtered TADs to find borders 
# and I set the boundaries limits based on changes in interactions. Note the boundaries are restricted to the inside of TADs.
# This version of the script splits the calculation of interactions and the boundaries definition into 2 separate tasks to makes it
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
for(c in c(22:1,"X")){  # All matrices can be loaded, but it is safer/more efficient to load only 1 matrix and send 1 job/chromosome.
  tmp <- gzfile(paste0("norm_hic_data/GM12878/chr",c,"_5kb_norm.txt.gz"),"rt")  # Files are unzipped before read
  matlist[[c]] <- c(readMM(tmp),c)  # Each entry in the list contains 2 object: the matrix and the chromosome name associated to the matrix
  close(tmp)  # Freeing connection for next matrix
}
TAD <- read.table("TAD/short_fullover_TAD.bed")  # Loading list of TADs
colnames(TAD) <- c("chr","start","end")  # More convenient column names (needed)

#==========================

#NOTE: To make the code faster, all 2D matrices have been converted to 1D vectors. This also makes the code more complex


#Defining functions

D2toD1 <- function(I,J,N){  # Simple function used to make the code more compact. Used to convert 2 dimension coordinates into 1D index
  # I: row number, J: column number, N: number of columns in the matrix
  return(((I-1)*N)+J)
}

vec_sub_square <- function(v,s,e,n,w){  
  # this function allows to get the values inside a square sub matrix (diamond) from a larger square matrix, however
  # all coordinates are in 1D (i.e. both objects are vectors).
  # v=vector,s=upper left corner of square, e=bottom right, n= number of cols in matrix,w=width of subsquare (diamond)
  sub_sq <- rep(0,times=w*w)  # submatrix that will hold the values
  cs <- 1  # counter used to iterate over indexes of the submatrix
  while(s<=e){  # as long as the current index is not higher than the last desired index (lower right corner of diamond)
    for(i in s:(s+(w-1))){  # Iterating over columns in each row
      sub_sq[cs] <-v[i]  # assigning each position in the diamond the corresponding value in the contact matrix.
      cs <- cs+1  # Incrementing of 1 position (i.e. incrementing from 1 column in the matrix)
    }
    s <- s+n  # Incrementing of n position (i.e. incrementing from 1 row in the matrix)
  }
  return(sub_sq) # Returning the diamond.
}


vec_diam_slide<-function(m,R=5000, D=100000){  #Vectorized version of the slider. MUCH FASTER!
  L <- length(m[1,])  # number of columns in the matrix
  M <- as.vector(t(m))  # transforming matrix into a vector
  diam <- rep(0,L) # preallocating space for diamond-summed data.
  diam[(D/R):(L-(D/R-1))] <- sapply(X = seq(from=D2toD1(D/R,D/R,L),to=D2toD1((L-(D/R-1)),(L-(D/R-1)),L),by=(L+1)),
                                    simplify = T, FUN= function(d){
                                      sum(vec_sub_square(v=M,s=(d-L*(D/R-1)),e=(d+(D/R-1)),n=L,w=(D/R)))})
  # Computing the sum of interactions in all diamonds of width w along the diagonal (i.e. sliding the diamond)
  # Note the diamond of width w cannot slide over the whole matrix. It needs to begin at w and stop at n-(w-1) 
  # otherwise it gets out of the matrix.
  return(diam)  # returning a vector containing the sum of interactions for each position along the diagonal
}

find_bound<-function(n,R=5000, D=100000){ # D is diamond size, R is resolution
  sub_TAD <- TAD[as.character(TAD$chr)==paste0("chr",n),]  # subsetting TAD dataframe by chromosome
  print(n)
  sub_TAD$Lbound.start = sub_TAD$Lbound.end = sub_TAD$Rbound.start =
    sub_TAD$Rbound.end = sub_TAD$maxint <- rep(0,length(sub_TAD$start))  # allocating space for max interaction observed in each TAD.
  diam <- scan(file = paste0("diam_sums/GM12878/chr",n,".txt"),what = numeric())  # Loading file containing vector of sums for all positions along the diagonal
  for(i in 1:length(sub_TAD$start)){  # Iterating over TADs
    start.ind <- pos_2_index(sub_TAD$start[i])  # transforming TAD start position to row in matrix. (TADs are already at 5kb res, so no need to approximate)
    end.ind <- pos_2_index(sub_TAD$end[i])  # Same for end position
    sub_TAD$maxint[i]<-max(diam[start.ind:end.ind])  # putting max interaction value for each TAD in dataframe.
    sl=el <- c(start.ind, diam[start.ind])  # Initiating left and right limits of left boundary

    sub_TAD$Lbound.start[i] <- index_2_pos(sl[1]) # writing the left boundary start coordinate into the output dataframe  
    while(el[2]<(diam[start.ind]+sub_TAD$maxint[i]/10)){  # End of left boundary -> sliding to the right
      el[1] <- el[1]+1; el[2] <-diam[el[1]]  # incrementing the end of the left boundary by 1 position
      if(el[1]==(length(diam)-(D/R))){break}}  # If the boundary is reaching the end of the part of the matrix over which the diamond could slide
    sub_TAD$Lbound.end[i] <- index_2_pos(el[1])  # writing the left boundary end coordinate into the output dataframe
    sr=er <-c(end.ind, diam[end.ind])  # Initiating left and right limits of right boundary
    
    while(sr[2]<(diam[end.ind]+sub_TAD$maxint[i]/10)){  # Start of right boundary -> sliding to the left
      sr[1] <- sr[1]-1; sr[2] <-diam[sr[1]]  # incrementing the start of the right boundary by 1 position
      if(sr[1]==pos_2_index(D)){break}} # If the boundary is reaching the end of the part of the matrix over which the diamond could slide
    sub_TAD$Rbound.start[i] <- index_2_pos(sr[1])  # writing the right boundary start coordinate into the output dataframe

    sub_TAD$Rbound.end[i] <- index_2_pos(er[1]) # writing the right boundary end coordinate into the output dataframe
  }
  return(sub_TAD)  # Returning a dataframe containing all TADs in the chromosome and their boundaries coordinates.
}

index_2_pos <- function(ind, resolution=5000) { return(resolution * (ind - 1)) }# Finding actual position in genome from bin number.
pos_2_index <- function(pos, resolution=5000) { return(1 + (pos / resolution))  }  # Finding index in vector from entry in matrix
#=========================

# 2 steps method:

# STEP 1: Only calculating the sums of interactions in each window and storing them: 
# If it still crashes, see alternative STEP 1 at the end to split the matrix.
for(n in c("9")){
  i <- matlist[[n]][[1]]
  diam<- vec_diam_slide(i)
  write.table(diam, file=paste0("diam_sums/GM12878/chr",matlist[[n]][[2]]),sep = "\t",quote=F,row.names = F,col.names = F)
  # Writing the vectors containing sums for each position into files so that matrices will not need to be loaded in the environment
  # when computing boundaries.
}

# STEP 2: loading summed interactions and defining boundaries from the vectors. No need to load contact matrices.
whole_bound <- data.frame()  # Dataframe used to store concatenated results
for(n in c(seq(1,22),"X")){
  bound<- find_bound(n)
  Lb <- bound[,c(1,9,8,4)]  # splitting dataframe into 2 separate dataframes: 1 for left boundaries and 1 for right boundaries
  Rb <- bound[,c(1,7,6,4)]
  colnames(Lb)=colnames(Rb) <- c("chr","start","end","ID")  # Making colnames match for the right and left dataframes
  whole_bound <- rbind(whole_bound,Lb)  # concatenating the 2 dataframes so that 1 boundary = 1 entry. 
  whole_bound <- rbind(whole_bound,Rb)
}
write.table(whole_bound,"TAD/hicboundaries/short_fullover/GM12878_HIC_boundaries.bed",col.names=F,row.names=F,quote=F,sep="\t")


##################################################
# Alternative STEP 1 in case the matrix is too large and needs to be splitted into overlapping parts:
# Used for chromosome 2, as its contact matrix was too heavy.
# Tested on chromosome 3, results were exactly identical to what is obtained with the regular method.

D <- 100000; R <- 5000; P <- D/R # The size of the sliding diamond and the resolution will determine how the 2 matrices should overlap.

for(n in c("3")){
  k <- dim(matlist[[n]][[1]])[1]  # number of rows in the matrix
  i <- matlist[[n]][[1]][1:(k/2 + P),1:(k/2 + P)]  # First subsetted matrix
  diam.i<- vec_diam_slide(i)  # computing sums in first subsetted matrix
  j <- matlist[[n]][[1]][(k/2 - P):k,(k/2 - P):k]  # second subsetted matrix
  diam.j<- vec_diam_slide(j)  # computing sums in second subsetted matrix
  diam <- append(diam.i[-((length(diam.i)-(P-1)):length(diam.i))],diam.j[-(1:(P+1))])  
  # Trimming the overlapping part of the output vectors for both subsets and merging the results.
  write.table(diam, file=paste0("diam_sums/GM12878/comp_chr",matlist[[n]][[2]]),sep = "\t",quote=F,row.names = F,col.names = F)
}
