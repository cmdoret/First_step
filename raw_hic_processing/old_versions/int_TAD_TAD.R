# This script computes the mean interaction in each TAD, given a matrix of contact and a list of TADs.
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
diam_slide<-function(m,R=5000,tad){  # L = vector containing each gene length; R = resolution, m = intrachromosomal matrix
  N <- length(m[[1]][1,])  # Storing number of cols
  M <- m[[1]]  # Transforming matrix into vector
  diam <- rep(0,times=length(tad[,1])) # preallocating space for diamond-summed data.
  c <- 1
  for(r in 1:length(tad$ID)){  # Iterating over TADs
    start <- pos_2_index(tad$start[r],R)
    end <- pos_2_index(tad$end[r],R)
    print(paste0("r=",r,"; N=",N,"; start=",start,"; end=",end))
    submat <- M[start:end,start:end]
    diam[c] <-  mean(submat)  # desired width of square (based on gene length) 
    # Storing normalized diamond sums in vector
    c <- c+1
  }
  return(cbind(tad,diam))
}

# Loading all matrices in a list (takes pretty long)
TAD_full <- read.table("TAD/short/short_fullover_TAD.bed")
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
  #linc <- linc_full[linc_full$chr==paste0("chr",c),] #only on the same chromosome as your matrix
  tmp_results <- diam_slide(matlist[[c]],tad=TAD) 
  print(tmp_results)
  results <- rbind(results,tmp_results)
}
rownames(results) <- NULL

#=========================
# Example call


setwd("/scratch/beegfs/monthly/mls_2016/cmatthey/first_step/staging_area/")
write.table(results, file = "TAD_contact/contact_linc_TAD_GM12878.txt",quote = F,col.names = T,row.names = F,sep = "\t")

