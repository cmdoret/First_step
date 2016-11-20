# This script computes the interaction in each genes given a list of genes and an interaction matrix
# Here, all interactions between genes and TAD are computed.
# Cyril Matthey-Doret
# 20.10.2016

############################
#ENTER CHROMOSOMES HERE
chrom <- c(seq(1,22),"X")
############################

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
  for(d in 1:length(L)){i[d]<-floor(pos_2_index(S[d]))}  #transforming start position into rounded vector indexes
  D<-floor(pos_2_index(L))  # Transforming length of genes into number of bins.
  int_tad <- Intervals(matrix(c(tad$start,tad$end),ncol=2))
  int_linc <- Intervals(matrix(c(S,S+L),ncol=2))
  rownames(int_tad) <- tad$ID 
  for(r in i){  # Iterating over genes
    over <-interval_overlap(int_tad,int_linc[c])
    gTAD <- tad[tad$ID==names(unlist(over[over!=0])),]
    if(length(gTAD[,1])>0){
      gTAD$start <- pos_2_index(gTAD$start)
      gTAD$end <- pos_2_index(gTAD$end)
      print(paste0("r=",r,"; N=",N,"; D[c]=", D[c]))
      diam[c] <-  sum(M[gTAD$start:gTAD$end,r:r+D[c]])*(L[c]/R)  # desired width of square (based on gene length) 
    } else{
      diam[c] <- 0
    }
    # Storing normalized diamond sums in vector
    c <- c+1
  }
  return(cbind(linc,diam))
}

# Loading all matrices in a list (takes pretty long)
linc_full <- read.table("jendata/eqtl.lincRNA.bed")  #your set of genes
colnames(linc_full)<-c("chr","start","end","gene","strand") #more convenient
TAD_full <- read.table("TAD/short/short_TADs.bed")
colnames(TAD_full)<-c("chr","start","end","ID") #more convenient
matlist <- list()
results <- data.frame()
for(c in c(chrom)){  # Loading matrices
  tmp <- gzfile(paste0("norm_hic_data/GM12878/chr",c,"_5kb_norm.txt.gz"),"rt")  # Files are unzipped before read
  matlist[[c]] <- c(readMM(tmp),c)
  close(tmp)  # Freeing connection for next matrix
}
for(c in c(chrom)){  # Calling function and storing results for all chromosomes
  TAD <- TAD_full[TAD_full$chr==paste0("chr",c),]
  linc <- linc_full[linc_full$chr==paste0("chr",c),] #only on the same chromosome as your matrix
  if(length(linc$gene)>0){
    tmp_results <- diam_slide(matlist[[c]], S=linc$start, L=as.vector(linc$end-linc$start),tad=TAD) 
    print(tmp_results)
    results <- rbind(results,tmp_results)
  }
}
rownames(results) <- NULL


#=========================
# Example call
# Exporting the list of TADs and the required functions to each node in the cluster.

setwd("/scratch/beegfs/monthly/mls_2016/cmatthey/first_step/staging_area/")
write.table(results, file = "TAD_contact/contact_linc_TAD_GM12878.txt",quote = F,col.names = T,row.names = F,sep = "\t")

