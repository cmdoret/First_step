# In this script, I try to define TAD boundaries using the insulation score as described in Smith et al, 2016.
# Cyril Matthey-Doret
# 20.10.2016
############################

# Loading data
setwd("/Users/cmatthe5/Documents/First_step/data/")
raw_hic <- read.table("raw_hic_chr1_5kbres/MAPQGE30/chr1_5kb.RAWobserved")
norm.vec<-read.table("raw_hic_chr1_5kbres/MAPQGE30/chr1_5kb.KRnorm")

# Normalizing matrix using the normalization vector provided with the Hi-C data.
# The vector was calculated using a matrix balancing algorithm.

# The matrix is under the parse format. Each line is a value, and 0 values are ignored.
# A row is "row", "col", "value"

# Normalizing the matrix: for each line: val/(norm.vec[row/res + 1] * norm.vec[col/res + 1])

normatrix <- function(i){
  r <- raw_hic[i,]
  nr <- r/(norm.vec[(r[1]/5000)+1] * norm.vec[(r[2]/5000)+1])
  return(nr)
}

apply(X = raw_hic, MARGIN = 1, FUN = normatrix, res=5000)
indexlist <- seq(1,length(raw_hic$V1))
normalized <-lapply(X = list(indexlist), FUN = normatrix)

df <- data.frame(a=c(1,2,3,4),b=c(10,20,30,40),c=c(100,200,300,400))
apply(X=df,MARGIN = c(1,2),FUN=normatrix)
