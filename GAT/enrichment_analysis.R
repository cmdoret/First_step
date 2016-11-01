# In this script, I analyse the results obtained from GAT.
# Cyril Matthey-Doret
# Tue Nov  1 20:07:39 2016 ------------------------------

# Loading data
setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
pc_epc.pr_bins <- read.table("GAT/out/gat_bins/results/all_pc/gat_W_allpc_S_epc_pr_A_short5bins.tsv",header=T)
pc_bins_epc.pr <- read.table("GAT/out/gat_bins/results/all_pc/gat_W_allpc_S_short5bins_A_epc_pr.tsv",header=T)
#Open all files in folder in a loop


# Use function on all files and process data into a single figure/table



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
    if(gat_args[[i]][7]=="A"){a=8;e=6} else{a=7;e=8}  # Accounting for different filename structures
    out <- data.frame(workspace=gat_args[[i]][3],  # Cenerating temporary dataframe to be concatenated.
                      segment=gat_args[[i]][5],
                      annotation=gat_args[[i]][a],
                      element = gat_args[[i]][e],  # Can be either pr (promoter) or prb (promoter + body)
                      track=as.character(tmp_df$track),
                      pval=tmp_df$pvalue,
                      qval=tmp_df$pvalue,
                      fold=tmp_df$fold)
    final <- rbind(final, out)  # Concatenates processed results table of each GAT test into the final dataframe.
  }
  return(final) 
}


gat_gat <- condense_gat("GAT/out/gat_bins/results/")
chipseq_gat <- condense_gat("GAT/out/gat_chipseq/results/")
