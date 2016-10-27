# This script "slides" a diamond of 100kb*100kb with a stepsize of 5kb on the normalized Hi-C contact matrix. 
# All the interactions in the diamond are summed at every step and stored in a dataframe that can then be used to define TAD boundaries.
# Cyril Matthey-Doret
# Wed Oct 26 20:48:51 2016 ------------------------------
######################################

# Loading data
#setwd("/Users/cmatthe5/Documents/First_step/data/")
#setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
setwd("/scratch/beegfs/monthly/mls_2016/cmatthey/first_step/data/")

stopifnot(require(Matrix))
# Loading all matrices in a list (takes pretty long)

c=1
tmp <- gzfile(paste0("norm_hic_data/chr",c,"_5kb_norm.txt.gz"),"rt")  # Files are unzipped before read
mat <- readMM(tmp)



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


vec_diam_slide<-function(m,R=5000, D=100000){  #Vectorized version of the slider. MUCH FASTER!
  L <- length(m[1,])
  M <- as.vector(t(m))
  diam <- rep(0,L) # preallocating space for diamond-summed data.
  diam[(D/R+1):(L-(D/R))] <- sapply(X = seq(from=(D/R+L*D/R+1),to=(L*L-(D/R+L*D/R)),by=(L+1)),
                                                simplify = T, FUN= function(d){
                                                  sum(vec_sub_square(v=M,s=(d-L*(D/R-1)),e=(d+(D/R-1)),n=L,w=(D/R)))})
  return(diam)
}

for_diam_slide<-function(m,R=5000, D=100000){  #Vectorized version of the slider. MUCH FASTER!
  L <- length(m[1,])
  M <- as.vector(t(m))
  diam <- rep(0,L) # preallocating space for diamond-summed data.
  c <- 1
  for(r in (D/R+1):(L-D/R)){
    i = seq(from=(D/R+L*D/R+1),to=(L*L-(D/R+L*D/R)),by=(L+1))  # Vector containing all diagonal positions
    diam[r] <- sum(vec_sub_square(v=M,s=(i[c]-L*(D/R-1)),e=(i[c]+D/R-1),n=L,w=(D/R))) 
    # For each position on the diagonal, summing all interactions inside a square with its lower left corner on the diagonal.
    c <- c+1
  }
  return(diam)
}


#===============================================================
setwd("/scratch/beegfs/monthly/mls_2016/cmatthey/first_step/staging_area")
write.table(data.frame(diam_sum=for_diam_slide(mat)))