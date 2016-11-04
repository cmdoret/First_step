
# script to normalise the raw HiC data
# NB: I gzipped the HiC data
#


#setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")  # Path on PC
setwd("/Users/cmatthe5/Documents/First_step/data/")  # Path on CHUV mac
for(c in c("9")){
  #chr <- paste0("/home/cyril/Documents/Master/sem_1/First_step/data/raw_hic_5kb_res/chr",c,"/MAPQGE30/chr",c)  # Path on PC
  chr <- paste0("/Users/cmatthe5/Documents/First_step/data/raw_hic_5kb_res/GM12878/chr",c,"/MAPQGE30/chr",c)  # Path on CHUV mac
  
  ## load data ## ====================================
  
  # load sparse matrix data
  library(data.table)                                    # allows the use of fread function. It's faster !
  raw_hic <- fread(
    paste0("gunzip < ",file.path(
      paste0(chr, "_5kb.RAWobserved.gz")),""),       # passing shell command to unzip file on the fly without affecting the original.
    sep="\t", header=FALSE)
  colnames(raw_hic) <- c("row_pos","col_pos", "value")   # more convenient
  raw_hic <- as.data.frame(raw_hic)                      # NB: raw_hic is a data.table, not a data.frame, so we'll convert it
  
  # load normalisation vector data
  gzfc <- gzfile(file.path(paste0(chr, "_5kb.SQRTVCnorm.gz")),"rt") # Loading .gz compressed vector
  norm.vec <- read.table(gzfc)  # reading vector
  close(gzfc)  # Closing connection with .gz file
  norm.vec <- norm.vec[,1] # make into vector object
  # NB: The normalization vector has one entry per row/col in the matrix
  
  
  
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
  # (if needed)
  # a sparse matrix will hardly take up any more memory than the data.frame
  # and MUCH less than a normal matrix IF it is, in fact, sparse (ie many zeros)
  library(Matrix)
  
  raw.hic.sm <- sparseMatrix(i = raw_hic$row_index, 
                               j = raw_hic$col_index, 
                               x = raw_hic$norm.value,
                               dims = c(length(norm.vec), length(norm.vec)),
                               dimnames=list(matrix.starnames, matrix.starnames))
  
  

  writeMM(obj = raw.hic.sm,file=paste0("norm_hic_data/GM12878","chr",c,"_5kb_norm.txt"))
  # Writing file as a matrix object (could be written as a dataframe)
  # Next step is to compute the sum of interactions for each bin.
}

# EOF