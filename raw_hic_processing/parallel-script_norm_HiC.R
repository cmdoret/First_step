#
# script to normalise the raw HiC data
# NB: I gzipped the HiC data
#


if (is.null(n.cores)) { 
  stopifnot(require(parallel)) 
  n.cores <- parallel:::detectCores() 
}
stopifnot(require(snow))

clus <- makeSOCKcluster( n.cores  )
on.exit( stopCluster(clus) )  
setwd("/Users/cmatthe5/Documents/First_step/data/")

norm<-function(c){
  
  
  ## load data ## ====================================
  
  # load sparse matrix data
  library(data.table)                                    # pour utiliser fread - ca lira plus vite!
  raw_hic <- fread(
    paste0("gunzip < ",file.path(
      paste0(chr, "_5kb.RAWobserved.gz")),""),       # on passe une commande shell pour dezipper on-the-fly sans que ca n'affecte le zip original
    sep="\t", header=FALSE)
  colnames(raw_hic) <- c("row_pos","col_pos", "value")   # pour reference plus facile
  raw_hic <- as.data.frame(raw_hic)                      # NB: raw_hic is a data.table, not a data.frame, so we'll convert it
  
  # load normalisation vector data
  gzfc <- gzfile(file.path(paste0(chr, "_5kb.KRnorm.gz")),"rt")
  norm.vec <- read.table(gzfc)
  close(gzfc)
  norm.vec <- norm.vec[,1] # make into vector
  # NB: il y a une entree par ligne/colonne de la matrice
  
  
  
  
  ## normalise values ## ==============================
  # functions to convert positions into indexes & vice-versa
  pos_2_index <- function(pos, resolution=5000) { return(1 + (pos / resolution))  }
  index_2_pos <- function(ind, resolution=5000) { return(resolution * (ind - 1))  }
  
  
  # list all the possible POSITIONS 
  # as suggested by the length of the normalisation vector
  # (since not all positions are actually "seen" in the raw_hic data.frame)
  matrix.starnames <- index_2_pos(1:length(norm.vec))
  
  
  # prepare normalisation factors & normalise
  raw_hic$row_index   <- pos_2_index(raw_hic$row_pos)
  raw_hic$col_index   <- pos_2_index(raw_hic$col_pos)
  raw_hic$norm.factor <- norm.vec[raw_hic$row_index] * norm.vec[raw_hic$col_index]
  raw_hic$norm.value  <- (raw_hic$value / raw_hic$norm.factor)
  
  
  
  
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
  
  
  #write.table(x=raw.hic.sm,file="norm_hic_data/chr1_5kb_norm.txt",quote=F,row.names = F,col.names = F,sep="\t")
  return(raw.hic.sm)
}

clusterCall( clus, function() { 
  require(Matrix) ; require(data.table)       # appel pour executer ce code sur chaque noeud. 
  # Utile au debut pour charger les librairies necessaires dans chaque instance de R.
} )
#clusterExport(clus, c("truc", "muche", ...???), envir=environment())

tmp <- clusterApplyLB( clus, c(1:22,"X"), norm)

# EOF