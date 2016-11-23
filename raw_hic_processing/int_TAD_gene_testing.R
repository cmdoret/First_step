# This script computes the interaction in each genes given a list of genes and an interaction matrix
# Here, all interactions between genes and TAD are computed.
# Cyril Matthey-Doret
# 20.10.2016


# Loading data
setwd("/Users/cmatthe5/Documents/First_step/data/")
setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
#setwd("/scratch/beegfs/monthly/mls_2016/cmatthey/first_step/data/")
#setwd("/home/cyril/Documents/First_step/data/")

stopifnot(require(Matrix))
stopifnot(require(intervals))
#==========================

#Defining functions

index_2_pos <- function(ind, resolution=5000) { return(resolution * (ind - 1)) }# Finding actual position in genome from bin number.
pos_2_index <- function(pos, resolution=5000) { return(1 + (pos / resolution))  }  # Finding index in vector from entry in matrix

# m=contact matrix; R=resolution; L=vector containing length of genes;S=vector containing start of genes
diam_slide<-function(m,R=5000,S , L, tad){  # L = vector containing each gene length; R = resolution, m = intrachromosomal matrix
  N <- length(m[[1]][1,])  # Storing number of cols
  M <- m[[1]]  # Transforming matrix into vector
  diam <- rep(0,times=length(L)) # preallocating space for diamond-summed data.
  c <- 1
  i <- rep(0,length(L)) #preallocating space for vector indexes of start sites
  for(d in 1:length(L)){i[d]<-floor(pos_2_index(S[d],R))}  #transforming start position into rounded vector indexes
  E <- ceiling(pos_2_index(S+L,R))  # storing indexes of end positions
  int_tad <- Intervals(matrix(c(tad$start,tad$end),ncol=2))
  int_linc <- Intervals(matrix(c(S,S+L),ncol=2))
  rownames(int_tad) <- tad$ID
  for(r in i){  # Iterating over genes
    over <-interval_overlap(int_tad,int_linc[c])
    gTAD <- tad[tad$ID==names(unlist(over[over!=0])),]
    if(length(gTAD[,1])>0){
      gTAD$start <- pos_2_index(gTAD$start,R)
      gTAD$end <- pos_2_index(gTAD$end,R)
      print(paste0("r=",r,"; N=",N,"; E[c]=", E[c]))
      submat <- M[gTAD$start:gTAD$end,r:E[c]]
      print(submat)
      first.last_avg <- rowMeans(submat[,c(1,length(submat[1,]))])
      submat[,c(1,length(submat[1,]))] <- first.last_avg
      diam[c] <-  mean(submat[,1:(length(submat[1,])-1)])  # desired width of square (based on gene length) 
    } else{
      diam[c] <- 0
    }
    # Storing normalized diamond sums in vector
    c <- c+1
  }
  return(cbind(linc,diam))
}


testmat <- list(matrix(rep(0,100),nrow=10),"1")
for(i in 1:10){testmat[[1]][i,i] <- 1}
for(i in 4:7){testmat[[1]][4:7,i]<-1}

# Loading all matrices in a list (takes pretty long)
linc <- data.frame(chr=c(1,1),start=c(14,16),end=c(16,22),gene=c("a","b"),strand=c("+","+"))
TAD <- data.frame(chr=1,start=15,end=30,ID="TAD")

diam_slide(testmat, S=linc$start, L=as.vector(linc$end-linc$start),tad=TAD,R = 5) 


#=========================
# Example call
# Exporting the list of TADs and the required functions to each node in the cluster.

setwd("/scratch/beegfs/monthly/mls_2016/cmatthey/first_step/staging_area/")
write.table(results, file = "TAD_contact/contact_linc_TAD_GM12878.txt",quote = F,col.names = T,row.names = F,sep = "\t")

