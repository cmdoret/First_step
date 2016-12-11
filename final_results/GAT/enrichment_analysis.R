# In this script, I condense into a table the results obtained from GAT for further visualisation.
# Cyril Matthey-Doret
# Tue Nov  1 20:07:39 2016 ------------------------------

# Loading data
setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
setwd("/Users/cmatthe5/Documents/First_step/data/")
setwd("/home/cyril/Documents/First_step/data/")







condense_gat <- function(path){
  # This function recursively concatenates the results of all gat 
  # output files in a folder. It takes a path containing all 
  # subdirectories as argument.
  filenames <- list.files(path,recursive = T)  # List all files recursively
  data_names <- gsub("[.]tsv", "", basename(filenames))  # Removes the extension and path from filenames and stores them in a vector
  gat_args <- strsplit(data_names,"_")  # Each filename is split by "_" into a vector.
  final <- data.frame()  # Initiating dataframe that will contain all concatenated results
  for(i in 1:length(filenames)){  # Iterating over files
    tmp_df <- read.table(file.path(path, filenames[i]),header=T)  # Generating temporary dataframe (raw gat output table)
    if(gat_args[[i]][8]=="A"){a=9;e=7} else{a=7;e=9}  # Accounting for different filename structures
    out <- data.frame(workspace=gat_args[[i]][3],  # Cenerating temporary dataframe to be concatenated.
                      segment=gat_args[[i]][5],
                      annotation=gat_args[[i]][a],
                      element = gat_args[[i]][e],  # Can be either pr (promoter) or prb (promoter + body)
                      track=as.character(tmp_df$track),
                      annottrack=as.character(tmp_df$annotation),
                      pval=tmp_df$pvalue,
                      qval=tmp_df$qvalue,
                      fold=tmp_df$fold)
    final <- rbind(final, out)  # Concatenates processed results table of each GAT test into the final dataframe.
  }
  return(final) 
}


gat_gat <- condense_gat("GAT/out/10k_samples/gat_hic_boundaries/out/segments_overlap/results")
chipseq_gat <- condense_gat("GAT/out/10k_samples/chipgat_fullover/segments_overlap/results/")
control_gat <- condense_gat("GAT/out/pos_control_gat/")
whole_GAT <- rbind(gat_gat,chipseq_gat)
write.table(whole_GAT,"GAT/out/whole_seg_10kgat_test_results.txt",quote=F,sep="\t",row.names=F)

