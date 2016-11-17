# In this script I look at the proportion of overlapping CTCF, RAD21 and SMC3
# Cyril Matthey-Doret
# Fri Nov 11 16:06:41 2016 ------------------------------

setwd("/Users/cmatthe5/Documents/First_step/data/")

# Loading bindingsites and overlaps
RAD21 <- read.table("chip_seq/GM12878_RAD21_peaks.bed")
CTCF <- read.table("chip_seq/GM12878_CTCF_peaks.bed")
SMC3 <- read.table("chip_seq/GM12878_SMC3_peaks.bed")
RAD21_CTCF <- read.table("chip_seq/inter_insulators_GM12878/inter_CTCF_RAD21.bed")
RAD21_SMC3 <- read.table("chip_seq/inter_insulators_GM12878/inter_RAD21_SMC3.bed")
SMC3_CTCF <- read.table("chip_seq/inter_insulators_GM12878/inter_CTCF_SMC3.bed")
all_inter <- read.table("chip_seq/inter_insulators_GM12878/inter_CTCF_RAD21_SMC3.bed")

# Building Venn diagram:
library(VennDiagram)

grid.newpage()
draw.triple.venn(euler.d = T,scaled = T, area1 = length(RAD21$V1), area2 = length(CTCF$V1), area3 = length(SMC3$V1), n12 = length(RAD21_CTCF$V1), 
                 n23 = length(SMC3_CTCF$V1), n13 = length(RAD21_SMC3$V1), 
                 n123 = length(all_inter$V1), category = c("RAD21", "CTCF", "SMC3"), lty = "blank", 
                 fill = c("skyblue", "yellow", "mediumorchid"))

# Making another one with proportional areas:
library("venneuler")
vd <- venneuler(c(RAD21 = length(RAD21$V1), CTCF = length(CTCF$V1), SMC3 = length(SMC3$V1), "RAD21&CTCF" = length(RAD21_CTCF$V1), 
                 "CTCF&SMC3" = length(SMC3_CTCF$V1), "RAD21&SMC3" = length(RAD21_SMC3$V1), 
                 "RAD21&CTCF&SMC3" = length(all_inter$V1)))
plot(vd)

vd <- venneuler(c(A=0.3, B=0.3, C=1.1, "A&B"=0.1, "A&C"=0.2, "B&C"=0.1 ,"A&B&C"=0.1))
plot(vd)
