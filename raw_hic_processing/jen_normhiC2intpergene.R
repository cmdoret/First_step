# This script computes the interaction in each genes given a list of genes and an interaction matrix
# Hi Jennifer, note you will need to adjust the paths and the code, depending on what chromosome you want.
# the variable containing the genes is called linc, but you can just use pc genes as well.
# Note that you need to only include genes in the chromosome of the matrix (I did that on line 27, you can change the chromosome from there)
# I hope I did not make any mistake ! Have fun :)
# Cyril Matthey-Doret
# 20.10.2016

############################

# Loading data
setwd("/Users/cmatthe5/Documents/First_step/data/")
#setwd("/scratch/beegfs/monthly/mls_2016/cmatthey/first_step/data/")
#setwd("/home/cyril/Documents/First_step/data/")

stopifnot(require(Matrix))
# Loading all matrices in a list (takes pretty long)
matlist <- list()
for(c in c("1")){
  tmp <- gzfile(paste0("norm_hic_data/chr",c,"_5kb_norm.txt.gz"),"rt")  # Files are unzipped before read
  matlist[[c]] <- c(readMM(tmp),c)
  close(tmp)  # Freeing connection for next matrix
}

linc <- read.table("jendata/eqtl.lincRNA.bed")  #your set of genes
colnames(linc)<-c("chr","start","end","gene","strand") #more convenient
linc <- linc[linc$chr=="chr1",] #only on the same chromosome as your matrix
#==========================


#Defining functions

vec_sub_square <- function(v,s,e,n,w){  
  # this function allows to get the values inside a square sub matrix in a 1D vector representing a larger square matrix
  # v=vector,s=upper left corner of square, e=bottom right, n= number of cols in matrix,w=width of subsquare
  sub_sq <- rep(0,times=w)
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

D2toD1 <- function(I,J,N){  #I=row,J=col,N=total number of columns or rows
  return(((I-1)*N)+J)
}

index_2_pos <- function(ind, resolution=5000) { return(resolution * (ind - 1)) }# Finding actual position in genome from bin number.
pos_2_index <- function(pos, resolution=5000) { return(1 + (pos / resolution))  }  # Finding index in vector from entry in matrix

# m=contact matrix; R=resolution; L=vector containing length of genes;S=vector containing start of genes
for_diam_slide<-function(m,R=5000,S , L){  #Vectorized version of the slider. MUCH FASTER! uses for loop instead of apply. Less greedy ?
  N <- length(m[[1]][1,])  # Storing number of cols
  M <- as.vector(t(m[[1]]))  # Transforming matrix into vector
  diam <- rep(0,times=length(L)) # preallocating space for diamond-summed data.
  c <- 1
  i <- rep(0,length(L)) #preallocating space for vector indexes of start sites
  for(d in 1:length(L)){i[d]<-D2toD1(floor(pos_2_index(S[d])),floor(pos_2_index(S[d])),N)}  #transforming start position into rounded vector indexes
  D<-ceiling(pos_2_index(L))  # Transforming length of genes into length of bins
  for(r in i){  # Iterating over genes
    diam[c] <- (L[c]/R)*sum(vec_sub_square(v=M,  # matrix
                                           s=(i[c]), # upper left corner of the square
                                           e=(i[c]+D[c]*N+D[c]), # lower right corner
                                           n=N,  # total n col of matrix
                                           w=(D[c])))  # desired width of square (based on gene length)
    # Storing normalized diamond sums in vector
    c <- c+1
  }
  return(cbind(linc,diam))
}
for_diam_slide(matlist[["1"]], S=linc$start, L=as.vector(linc$end-linc$start))



#=========================
# Example call
# Exporting the list of TADs and the required functions to each node in the cluster.

full <- for_diam_slide(matlist[["1"]], S=linc$start, L=as.vector(linc$end-linc$start)) # The diamond size depends on the gene length! it is not fixed
setwd("/scratch/beegfs/monthly/mls_2016/cmatthey/first_step/staging_area/")
write.table(full, file = "chr21_linc.txt",quote = F,col.names = T,row.names = F,sep = "\t")

