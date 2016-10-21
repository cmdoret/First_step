# In this script, I normalize the hicdata matrix.
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
write.table(x = normalized,file = "writetable.raw_observed_norm_chr1_5kb.txt",sep = "\t",quote = F,row.names = F,col.names = F)



