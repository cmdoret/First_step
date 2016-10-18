# The purpose of this script is to compute the median length of overlap with RNAs for each TAD bin.
# The bins have been defined with three different width (1%, 5% and 10% of their TAD length).
# Cyril Matthey-Doret
# 17.10.2016
###################################

#Loading data:
#setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
setwd("/home/cyril/Documents/First_step/data/")
over.len_pc_bins<-read.table("TAD/merged/overlaps/length_bin-pcgenes_overlap10.bed")
over.len_lincRNA_bins<-read.table("TAD/merged/overlaps/length_bin-lincRNA_overlap10.bed")
colnames(over.len_pc_bins)=colnames(over.len_lincRNA_bins) <- c("bin.chr","bin.start","bin.end","bin.ID","bin.pos",
                                                                "RNA.chr","RNA.start","RNA.end","RNA.ID","RNA.strand","overlap.length")
#=================================

# Normalizing per bin length
over.len_pc_bins <- data.frame(over.len_pc_bins, norm.over=over.len_pc_bins$overlap.length/(over.len_pc_bins$RNA.end - over.len_pc_bins$RNA.start))
over.len_lincRNA_bins <- data.frame(over.len_lincRNA_bins,norm.over=over.len_lincRNA_bins$overlap.length/(over.len_lincRNA_bins$RNA.end - over.len_lincRNA_bins$RNA.start))
      
medover.lincRNA <-by(data = over.len_lincRNA_bins$norm.over[over.len_lincRNA_bins$overlap.length>0],
   INDICES = over.len_lincRNA_bins$bin.pos[over.len_lincRNA_bins$overlap.length>0],
   FUN = median,simplify = F)

medover.pcgene <-by(data = over.len_pc_bins$norm.over[over.len_pc_bins$overlap.length>0],
                        INDICES = over.len_pc_bins$bin.pos[over.len_pc_bins$overlap.length>0],
                        FUN = median,simplify = F)

#===============================
# Visualize results:
thr <- 1/length(medover.pcgene) + 0.6/length(medover.pcgene)
medover_bin_RNA <- data.frame(pos=factor(names(unlist(medover.pcgene)),
                                       c(paste0("L",seq(from=(0.3/thr),to=1,by=-1)),
                                         seq(from=thr*100,to=100,by=thr*100),
                                         paste0("R",seq(from=1,to=(0.3/thr),by=1))),ordered = T),
                            median.pc.over_med=unname(unlist(medover.pcgene)),
                            median.lincRNA.over_med=unname(unlist(medover.lincRNA)))

par(mfrow=c(2,1))
plot(medover_bin_RNA$pos,medover_bin_RNA$median.lincRNA.over_med)
plot(medover_bin_RNA$pos,medover_bin_RNA$median.pc.over_med)

write.table(medover_bin_RNA,file = "TAD/merged/median_overlap_bin10-RNA.txt",sep="\t",quote = F,col.names = F,row.names = F)


