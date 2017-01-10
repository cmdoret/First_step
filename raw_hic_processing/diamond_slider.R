# This script "slides" a diamond of 100kb*100kb with a stepsize of 5kb on the normalized Hi-C contact matrix. 
# All the interactions in the diamond are summed at every step and stored in a dataframe that can then be used to define TAD boundaries.
# Cyril Matthey-Doret
# Wed Oct 26 20:48:51 2016 ------------------------------
######################################

# Loading data
#setwd("/Users/cmatthe5/Documents/First_step/data/")
#setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
setwd("/home/cyril/Documents/First_step/data/")
if (!exists("n.cores")) { 
  stopifnot(require(parallel)) 
  n.cores <- parallel:::detectCores()-1 # Value depends on number of cores in the CPU
}
stopifnot(require(snow))
stopifnot(require(Matrix))
# Loading all matrices in a list (takes pretty long)
matlist <- list()
for(c in c(22)){
  tmp <- gzfile(paste0("norm_hic_data/GM12878//chr",c,"_5kb_norm.txt.gz"),"rt")  # Files are unzipped before read
  matlist[[c]] <- c(readMM(tmp),c)
  close(tmp)  # Freeing connection for next matrix
}
# Converting matrix coordinates to vector index
D2toD1 <- function(I,J,N){
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


diam_slide<-function(m,R=5000, D=100000){
  M <- m
  diam <- rep(0,length(M[1,])) # preallocating space for diamond-summed data.
  diam[(D/R):(length(M[1,])-(D/R-1))] <- sapply(X = seq(D/R,(length(M[1,])-(D/R-1))),
                                                simplify = T, FUN= function(d){
                                                  sum(M[d:(d-(D/R-1)),d:(d+(D/R-1))])})
  return(diam)
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

for_diam_slide<-function(m,R=5000, D=100000){  #Vectorized version of the slider. MUCH FASTER!
  L <- length(m[1,])
  M <- as.vector(t(m))
  diam <- rep(0,L) # preallocating space for diamond-summed data.
  c <- 1
  i = seq(from=D2toD1(D/R,D/R,L),to=D2toD1((L-(D/R-1)),(L-(D/R-1)),L),by=(L+1))
  for(r in (D/R):(L-(D/R-1))){
    #print(i[c])
    diam[r] <- sum(vec_sub_square(v=M,s=(i[c]-L*(D/R-1)),e=(i[c]+D/R),n=L,w=(D/R)))
    c <- c+1
  }
  return(diam)
}

#===============================================================
# Start cluster:
clusterCall( clus, function() { 
  require(Matrix)       # executed on each node. 
  # Loading required libraries in each instance of R
} )
clusterExport(clus, c("vec_sub_square"), envir=environment())
tmp <- clusterApplyLB( clus, matlist, vec_diam_slide)

#===============================================================
#Benchmark
time_rec <- data.frame(size= seq(from=40,to=10000,by=10),time_mat=rep(0),time_vec=rep(0),res_mat=rep(0),res_vec=rep(0),
                       res_for=rep(0),time_for=rep(0))
for(i in seq(from=40,to=10000,by=10)){
  print(i)
  tpm <- proc.time()
  time_rec$res_mat[time_rec$size==i] <- sum(diam_slide(matlist[[22]][[1]][4000:(4000+i),4000:(4000+i)]))
  time_rec$time_mat[time_rec$size==i] <-(proc.time()-tpm)[3]
  tpm <- proc.time()
  time_rec$res_vec[time_rec$size==i] <- sum(vec_diam_slide(matlist[[22]][[1]][4000:(4000+i),4000:(4000+i)]))
  time_rec$time_vec[time_rec$size==i] <-(proc.time()-tpm)[3]
  tpm <- proc.time()
  time_rec$res_for[time_rec$size==i] <- sum(for_diam_slide(matlist[[22]][[1]][4000:(4000+i),4000:(4000+i)]))
  time_rec$time_for[time_rec$size==i] <-(proc.time()-tpm)[3]
}
plot(time_rec$size,time_rec$res_mat,xlab="Matrix size",ylab="Sum of results",
     main="speed comparison",col="red")
points(time_rec$size,time_rec$res_vec,col="blue")
points(time_rec$size,time_rec$res_for,col="green")
plot(time_rec$size,time_rec$time_mat,xlab="Matrix size"
     ,ylab="Time [s]",main="Diamond sliders benchmarking",col="red")
points(time_rec$size,time_rec$time_vec,col="blue")
points(time_rec$size,time_rec$time_for,col="green")
legend(fill = c("red","blue","green"), "topleft",legend=
         c("matrix structure: \"apply\"", "vectorized matrix: \"apply\"", "vectorized matrix: \"for loop\""))




a <- vec_diam_slide(matlist[[22]][[1]][4000:(4000+i),4000:(4000+i)])
b <- diam_slide(matlist[[22]][[1]][4000:(4000+i),4000:(4000+i)])
  
plot(a,b)

testmat <- list(matrix(rep(0,100),nrow=10),"1")
for(i in 1:10){testmat[[1]][i,i] <- 1}
for(i in 5:7){testmat[[1]][4,i]<-1}
for(i in 6:7){testmat[[1]][5,i]<-1}
testmat[[1]][6,7]<-1
