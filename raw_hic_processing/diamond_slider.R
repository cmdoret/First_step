# This script "slides" a diamond of 100kb*100kb with a stepsize of 5kb on the normalized Hi-C contact matrix. 
# All the interactions in the diamond are summed at every step and stored in a dataframe that can then be used to define TAD boundaries.
# Cyril Matthey-Doret
# Wed Oct 26 20:48:51 2016 ------------------------------
######################################

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


# Converting matrix coordinates to vector index
D2toD1 <- function(I,J,N){
  return(((I-1)*N)+J)
}

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

diam_slide<-function(m,R=5000, D=100000){
  M <- m
  diam <- rep(0,length(M[1,])) # preallocating space for diamond-summed data.
  diam[(D/R+1):(length(M[1,])-(D/R))] <- sapply(X = seq(D/R+1,(length(M[1,])-(D/R))),
                                                simplify = T, FUN= function(d){
                                                  sum(M[d:(d-(D/R+1)),d:(d+(D/R-1))])})
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
#Benchmark
time_rec <- data.frame(size= seq(from=40,to=1000,by=10),time_mat=rep(0),time_vec=rep(0),res_mat=rep(0),res_vec=rep(0))
for(i in seq(from=40,to=1000,by=10)){
  tpm <- proc.time()
  time_rec$res_mat[time_rec$size==i] <- sum(diam_slide(matlist[[22]][[1]][4000:(4000+i),4000:(4000+i)]))
  time_rec$time_mat[time_rec$size==i] <-(proc.time()-tpm)[3]
  tpm <- proc.time()
  time_rec$res_vec[time_rec$size==i] <- sum(vec_diam_slide(matlist[[22]][[1]][4000:(4000+i),4000:(4000+i)]))
  time_rec$time_vec[time_rec$size==i] <-(proc.time()-tpm)[3]
}
plot(time_rec$size,time_rec$res_mat,xlab="Number of iterations",ylab="Sum of results",
     main="Number of iterations \n(matrix increase of 10x10 per iteration)",col="red")
points(time_rec$size,time_rec$res_vec,col="blue")
plot(time_rec$size,time_rec$time_mat,xlab="Number of iterations \n(matrix increase of 10x10 per iteration)"
     ,ylab="Time [s]",main="Number of iterations on a matrix",col="red")
points(time_rec$size,time_rec$time_vec,col="blue")

# ISSUE: VEC_DIAM_SLIDE RETURNS A VECTOR OF LENGTH 21 instead of 21^2
