# In this script, I analyse the results obtained from GAT.
# Cyril Matthey-Doret
# Tue Nov  1 20:07:39 2016 ------------------------------

# Loading data
setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
setwd("/Users/cmatthe5/Documents/First_step/data/")
setwd("/home/cyril/Documents/First_step/data/")
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


gat_gat <- condense_gat("GAT/out/10k_samples/bingat_fullover/segments_overlap/results")
chipseq_gat <- condense_gat("GAT/out/10k_samples/chipgat_fullover/segments_overlap/results/")
control_gat <- condense_gat("GAT/out/pos_control_gat/")
whole_GAT <- rbind(gat_gat,chipseq_gat)
write.table(whole_GAT,"GAT/out/whole_seg_10kgat_test_results.txt",quote=F,sep="\t",row.names=F)

# Data analysis
whole_gat <- read.table("GAT/out/whole_seg_10kgat_results.txt",header=T)

bins2gene<-droplevels(whole_gat[whole_gat$segment=="short5bins",])
bins2gene$track <- as.character(bins2gene$track)
bins2gene$track <- factor(bins2gene$track,levels = c(1:20))

gene2bins<-droplevels(whole_gat[whole_gat$annotation=="short5bins",])
gene2bins$annottrack <- as.character(gene2bins$annottrack)
gene2bins$annottrack <- factor(gene2bins$annottrack,levels = c(1:20))

bs2gene<-droplevels(whole_gat[whole_gat$segment=="bindingsite",])
gene2bs<-droplevels(whole_gat[whole_gat$annotation=="bindingsite",])


library(ggplot2)
elements_names <- c(
  `pr` = "promoter",
  `prb` = "promoter + gene body",
  `elinc` = "elincRNA",
  `nelinc` = "non-elincRNA",
  `epc` = "epc gene",
  `nepc` = "non-epc gene"
)
b2g <- ggplot(data=bins2gene[bins2gene$workspace=="wholegenome",],aes(x=track,y=fold,fill=log10(qval)))+
  scale_fill_continuous(low="#66bb66",high="#bb6666",guide = "colourbar",name= "Log10 q-value")+
  geom_bar(stat="identity")+theme_bw()+
  facet_grid(element~annotation,labeller = as_labeller(elements_names))+
  ggtitle("Enrichment of TAD bins in lincRNAs and pcgenes")+
  xlab("TAD bins") + ylab("Fold enrichment")

g2b <- ggplot(data=gene2bins[gene2bins$workspace=="wholegenome",],aes(x=annottrack,y=fold,fill=log10(qval)))+
  scale_fill_continuous(low="#66bb66",high="#bb6666",guide = "colourbar",name= "Log10 q-value")+
  geom_bar(stat="identity")+theme_bw()+
  facet_grid(element~segment,labeller = as_labeller(elements_names))+
  ggtitle("Enrichment of lincRNAs and pc genes in TAD bins")+
  xlab("TAD bins") + ylab("Fold enrichment")

bs2g <- ggplot(data=bs2gene,aes(x=workspace,y=fold,fill=log10(qval)))+
  scale_fill_continuous(low="#66bb66",high="#bb6666",guide = "colourbar",name= "Log10 q-value")+
  geom_bar(stat="identity")+theme_bw()+
  facet_grid(element~annotation,labeller = as_labeller(elements_names))+
  ggtitle("Enrichment of CTCF/cohesin binding sites in lincRNAs and pcgenes")+
  xlab("Workspace") + ylab("Fold enrichment")+theme(axis.text.x=element_text(angle = -90, hjust = 0))

g2bs <- ggplot(data=gene2bs,aes(x=workspace,y=fold,fill=log10(qval)))+
  scale_fill_continuous(low="#66bb66",high="#bb6666",guide = "colourbar",name= "Log10 q-value")+
  geom_bar(stat="identity")+theme_bw()+
  facet_grid(element~segment,labeller = as_labeller(elements_names))+
  ggtitle("Enrichment of lincRNAs and pc genes in CTCF/cohesin binding sites")+
  xlab("Workspace") + ylab("Fold enrichment")+theme(axis.text.x=element_text(angle = -90, hjust = 0))