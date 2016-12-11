# This script computes the mean interaction per 5kb in each TAD given a list of genes, a list of TADs and a (normalized) interaction matrix.
# If a gene overlaps 2 different TADs, the interactions will be computed independently for each TAD and the gene will appear twice in the list.
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
# Note: Can only be run for 1 chromosome at a time!
diam_slide<-function(m,R=5000,S , L, tad){  # L = vector containing each gene length; R = resolution, m = intrachromosomal matrix
  N <- length(m[[1]][1,])  # Storing number of cols
  M <- m[[1]]  # Transforming matrix into vector
  diam <- rep(0,times=length(L)) # preallocating space for diamond-summed data.
  c <- 1  # counter used to iterate over genes
  i <- rep(0,length(L)) #preallocating space for vector indexes of start sites
  for(d in 1:length(L)){i[d]<-floor(pos_2_index(S[d],R))}  #transforming start position into rounded vector indexes
  E <- floor(pos_2_index(S+L,R))  # storing indexes of end positions
  int_tad <- Intervals(matrix(c(tad$start,tad$end),ncol=2))  # Transforming all TADs in the chromosome into intervals
  int_linc <- Intervals(matrix(c(S,S+L),ncol=2))  # Transforming all genes in the chromosome into intervals
  rownames(int_tad) <- tad$ID   # Giving TAD IDs as interval names
  for(r in i){  # Iterating over genes
    over <-interval_overlap(int_tad,int_linc[c])  # Checking for overlap between gene c and all TADs in the chromosome
    gTAD <- tad[tad$ID==names(unlist(over[over!=0])),]  # Storing coordinates of all TADs overlapping gene c
    print(length(gTAD[,1]))  # printing TAD coordinates (indication of progress, can be removed)
    if(length(gTAD[,1])>0){  # If the gene overlap at least 1 TAD
      gTAD$start <- pos_2_index(gTAD$start,R)  # Convert the start positions of all TADs overlapping gene c to matrix coordinates
      gTAD$end <- pos_2_index(gTAD$end,R) # Convert the end positions of all TADs overlapping gene c to matrix coordinates
      print(paste0("gene.start=",r,"; gene.end=", E[c],"; Matrix dimensions=",N)) # Debugging purpose, can be commented
      print(paste0("gTADstart=", gTAD$start[1], ": gTAD$end=", gTAD$end[1]))  # Debugging/Progress purpose
      submat <- M[gTAD$start[1]:gTAD$end[1],gTAD$start[1]:gTAD$end[1]]  # Generating a square sub-matrix spanning from the start to the end of the first TAD
      #print(submat)
      diam[c] <-  mean(submat)  # Storing the mean of all values in that submatrix as a value in the output vector
      if(length(gTAD[,1])>1){   # If more than 1 TAD overlap the gene c
        print(paste("DOUBLE Overlap", gTAD, SEP=": "))  # Print an information message stating there are multiple overlaps, can be removed
        for(t in 2:length(gTAD$start)){  # Iterating from second TAD to last
          print(paste0("r=",r,"; N=",N,"; E[c]=", E[c])) # Information message, can be removed
          print(paste0("gTAD$start=", gTAD$start[t], ": gTAD$end=", gTAD$end[t])) # Same here
          submat <- M[gTAD$start[t]:gTAD$end[t],gTAD$start[t]:gTAD$end[t]] # Generating a square submatrix for each overlapping TAD beyond the first
          linc <- rbind(linc,linc[c,]) # appending a duplicate entry for gene c to the input list of gene
          diam <-  append(diam,mean(submat))  # appending a new value at the end of the output vector that will match duplicate entry in input list
        }
      }
    } else{  # If the gene does not overlap any TAD
      diam[c] <- 0  # Do not count interactions (return 0)
    }
    # Storing normalized diamond sums in vector
    c <- c+1  # Incrementation to next gene
  }
  return(cbind(linc,diam))  # Binding mean interaction per TAD and gene list together and returning the dataframe
}

# Loading all matrices in a list (takes pretty long)
linc_full <- read.table("linc_RNA/LCL.expressed.lincRNA.bed")  # set of input genes
colnames(linc_full)<-c("chr","start","end","gene","strand")  # more convenient colnames (needed, the function will use these)
TAD_full <- read.table("TAD/short/short_fullover_TAD.bed")  # set of input TADs
colnames(TAD_full)<-c("chr","start","end","ID") #more convenient colnames
matlist <- list()  # Matrices are stored in a list
results <- data.frame()  # Will store the concatenated results for all chromosomes
for(c in c(chrom)){  # Loading matrices
  tmp <- gzfile(paste0("norm_hic_data/GM12878/chr",c,"_5kb_norm.txt.gz"),"rt")  # Files are unzipped before read
  matlist[[c]] <- c(readMM(tmp),c)  # Reading sparse matrices
  close(tmp)  # Freeing connection for next matrix
}
for(c in c(chrom)){  # Calling function and storing results for all chromosomes
  TAD <- TAD_full[TAD_full$chr==paste0("chr",c),] # subsetting TAD list for each chromosome
  linc <- linc_full[linc_full$chr==paste0("chr",c),] # subsetting gene list for each chromosome
  rownames(linc) <- NULL  # Resetting rownames to avoid problems
  if(length(linc$gene)>0){  # If there is at least 1 gene in the chromosome
    tmp_results <- diam_slide(matlist[[c]], S=linc$start, L=as.vector(linc$end-linc$start),tad=TAD) # Call function on each chromosome
    #print(tmp_results)
    results <- rbind(results,tmp_results)  # Concatenate results for each chromosome
  }
}
rownames(results) <- NULL  # Reset rownames


#=========================

setwd("/scratch/beegfs/monthly/mls_2016/cmatthey/first_step/staging_area/")
write.table(results, file = "TAD_contact/contact_linc_TAD_GM12878.txt",quote = F,col.names = T,row.names = F,sep = "\t")

