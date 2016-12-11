# This script is used to normalize the raw Hi-C data. The normalization follows the procedure described 
# in the README of Rao et al. 2014 supplementary data.(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525)
# Note: All raw matrices were gzipped beforehand.
# Cyril Matthey-Doret
# Sun Dec 11 20:13:19 2016 ------------------------------



#setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")  # Path on PC
setwd("/Users/cmatthe5/Documents/First_step/data/")  # Path on CHUV mac
for(c in c(seq(1:22),"X")){
  #chr <- paste0("/home/cyril/Documents/Master/sem_1/First_step/data/raw_hic_5kb_res/chr",c,"/MAPQGE30/chr",c)  # Path on PC
  chr <- paste0("/Users/cmatthe5/Documents/First_step/data/raw_hic_5kb_res/GM12878/chr",c,"/MAPQGE30/chr",c)  # Path on CHUV mac
  
  ## load data ## ====================================
  
  # load sparse matrix data
  library(data.table)                                    # allows the use of fread function
  raw_hic <- fread(
    paste0("gunzip < ",file.path(
      paste0(chr, "_5kb.RAWobserved.gz")),""),       # passing shell command to unzip file on the fly without affecting the original.
    sep="\t", header=FALSE)
  colnames(raw_hic) <- c("row_pos","col_pos", "value")   # more explicit column names
  raw_hic <- as.data.frame(raw_hic)                      # Converting raw_hic from data.table object to data.frame
  
  # load normalization vector data
  gzfc <- gzfile(file.path(paste0(chr, "_5kb.SQRTVCnorm.gz")),"rt") # Loading .gz compressed vector
  norm.vec <- read.table(gzfc)  # reading vector
  close(gzfc)  # Closing connection with .gz file
  norm.vec <- norm.vec[,1] # make into vector object
  # Note: The normalization vector has one entry per row/col in the matrix
  
  
  
  ## normalise values ## ==============================
  # functions to convert positions into indexes & vice-versa
  pos_2_index <- function(pos, resolution=5000) { return(1 + (pos / resolution))  }  # Finding index in vector from entry in matrix
  index_2_pos <- function(ind, resolution=5000) { return(resolution * (ind - 1))  }  # Finding entry in matrix corresponding to index in vector
  
  
  # list all the possible POSITIONS 
  # as suggested by the length of the normalisation vector
  # (since not all positions are actually "seen" in the raw_hic data.frame)
  matrix.starnames <- index_2_pos(1:length(norm.vec))
  
  
  # prepare normalisation factors & normalise
  raw_hic$row_index   <- pos_2_index(raw_hic$row_pos)
  raw_hic$col_index   <- pos_2_index(raw_hic$col_pos)
  raw_hic$norm.factor <- norm.vec[raw_hic$row_index] * norm.vec[raw_hic$col_index]
  raw_hic$norm.value  <- (raw_hic$value / raw_hic$norm.factor)
  raw_hic <- raw_hic[!is.na(raw_hic$norm.value),]
  
  
  
  ## prepare sparse matrix ## ========================
  # a sparse matrix will take much less space than a full matrix since Hi-C data is very sparse.
  library(Matrix)
  
  raw.hic.sm <- sparseMatrix(i = raw_hic$row_index, 
                               j = raw_hic$col_index, 
                               x = raw_hic$norm.value,
                               dims = c(length(norm.vec), length(norm.vec)),
                               dimnames=list(matrix.starnames, matrix.starnames))
  
  

  writeMM(obj = raw.hic.sm,file=paste0("norm_hic_data/GM12878","chr",c,"_5kb_norm.txt"))
}

# EOF