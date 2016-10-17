# The purpose of this script is to compute the median length of overlap with RNAs for each TAD bin.
# The bins have been defined with three different width (1%, 5% and 10% of their TAD length).
# Cyril Matthey-Doret
# 17.10.2016
###################################

#Loading data:
setwd("/home/cyril/Documents/Master/sem_1/First_step/data/")
over.len_pc_bins<-read.table("TAD/merged/overlaps/length_bin-pcgenes_overlap5.bed")
over.len_lincRNA_bins<-read.table("TAD/merged/overlaps/length_bin-lincRNA_overlap5.bed")
colnames(over.len_pc_bins)=colnames(over.len_lincRNA_bins) <- c("bin.chr","bin.start","bin.end","bin.ID","bin.pos",
                                                                "RNA.chr","RNA.start","RNA.end","RNA.ID","RNA.strand","overlap.length")
#=================================



# Should I include all overlaps in the statistics, or only those >0 ?
# Let's try both

# Without 0's:

medover.lincRNA_wo <-by(data = over.len_lincRNA_bins$overlap.length[over.len_lincRNA_bins$overlap.length>0],
   INDICES = over.len_lincRNA_bins$bin.pos[over.len_lincRNA_bins$overlap.length>0],
   FUN = median,simplify = F)

medover.pcgene_wo <-by(data = over.len_pc_bins$overlap.length[over.len_pc_bins$overlap.length>0],
                        INDICES = over.len_pc_bins$bin.pos[over.len_pc_bins$overlap.length>0],
                        FUN = median,simplify = F)

#with 0's:


medover.lincRNA_w <-by(data = over.len_lincRNA_bins$overlap.length,
                        INDICES = over.len_lincRNA_bins$bin.pos,
                        FUN = median,simplify = F)

medover.pcgene_w <-by(data = over.len_pc_bins$overlap.length,
                       INDICES = over.len_pc_bins$bin.pos,
                       FUN = median,simplify = F)

#===============================
# Visualize results:
thr <- 1/length(medover.pcgene_w) + 0.6/length(medover.pcgene_w)
medover_bin_RNA <- data.frame(pos=factor(names(unlist(medover.pcgene_w)),
                                       c(paste0("L",seq(from=(0.3/thr),to=1,by=-1)),
                                         seq(from=thr*100,to=100,by=thr*100),
                                         paste0("R",seq(from=1,to=(0.3/thr),by=1))),ordered = T),
                            median.pc.w=unname(unlist(medover.pcgene_w)),
                            median.pc.wo=unname(unlist(medover.pcgene_wo)),
                            median.lincRNA.w=unname(unlist(medover.lincRNA_w)),
                            median.lincRNA.wo=unname(unlist(medover.lincRNA_wo)))

plot(medover_bin_RNA$pos,medover_bin_RNA$median.pc.wo)
plot(medover_bin_RNA$pos,medover_bin_RNA$median.lincRNA.wo)

