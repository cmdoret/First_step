# In this script, I try to define TAD boundaries using the insulation score as described in Smith et al, 2016.
# Cyril Matthey-Doret
# 20.10.2016
#BSUB -o raw_observed_norm_chr1_5kb.txt
#BSUB -e hicnorm_errors.txt
#BSUB -J hic_normalization

############################

# Loading data
setwd("/Users/cmatthe5/Documents/First_step/data/")
raw_hic <- read.table("raw_hic_chr1_5kbres/MAPQGE30/chr1_5kb.RAWobserved")
norm.vec<-read.table("raw_hic_chr1_5kbres/MAPQGE30/chr1_5kb.KRnorm")
norm.vec <- as.vector(c(norm.vec))
# Normalizing matrix using the normalization vector provided with the Hi-C data.
# The vector was calculated using a matrix balancing algorithm.

# The matrix is under the parse format. Each line is a value, and 0 values are ignored.
# A row is "row", "col", "value"

# Normalizing the matrix: for each line: val/(norm.vec[row/res + 1] * norm.vec[col/res + 1])

normatrix <- function(i){
  r <- raw_hic[i,]
  nr <- r[3]/(norm.vec[[1]][as.numeric((r[1]/5000)+1)] * norm.vec[[1]][as.numeric((r[2]/5000)+1)])
  return(nr)
}

indexlist <- seq(1,length(raw_hic$V1))
normalized <-lapply(X = indexlist, FUN = normatrix)
normalized <-(unname(unlist(normalized)))
normalized <-cbind(raw_hic[,c(1,2)],norm=normalized)
normalized
write.table(x = normalized,file = "writetable.raw_observed_norm_chr1_5kb.txt",sep = "\t",quote = F,row.names = F,col.names = F)
#==================
# Tests
df <- data.frame(a=c(1,2,3,4),b=c(10,20,30,40),c=c(100,200,300,400))
apply(X=df,MARGIN = c(1,2),FUN=normatrix)

testfunc <- function(i){
  print(str(as.numeric(i)))
  return(i)
}

testlist <- seq(1,30)
testres<-lapply(X = testlist, FUN=testfunc)


