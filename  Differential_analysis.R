library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(GenomicRanges)
library(chromVAR) ## For FRiP analysis and differential analysis
library(DESeq2) ## For differential analysis section
library(ggpubr) ## For customizing figures
library(corrplot) ## For correlation plot 

# Set work place
setwd("/Users/lixia/Data/data/Cut_tag")
## Path to the project and histone list
projPath = "/Users/lixia/Data/data/Cut_tag"
histL = c("K27me3", "K4me3")
repL = c("rep1", "rep2")

### ====================================  Create a master peak list merging all the peaks called for each sample========================= ###
mPeak = GRanges()
## overlap with bam file to get count
for(hist in histL){
  for(rep in repL){
    peakRes = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_", rep, "_seacr_control.peaks.stringent.bed"), header = FALSE, fill = TRUE)
    mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
  }
}
masterPeak = reduce(mPeak)
head(masterPeak)

### ==================================== Get the fragment counts for each peak in the master peak list =========================== ###
bamDir = paste0(projPath, "/alignment/bam")
countMat = matrix(NA, length(masterPeak), length(histL)*length(repL))
## overlap with bam file to get count
i = 1
for(hist in histL){
  for(rep in repL){
    
    bamFile = paste0(bamDir, "/", hist, "_", rep, "_bowtie2.mapped.bam")
    fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts)[,1]
    i = i + 1
  }
}
colnames(countMat) = paste(rep(histL, each = length(repL)), rep(repL, length(histL)), sep = "_")

# Add chr coordinates  into rownames
rowRanges_names = paste0(seqnames(masterPeak), ":", start(masterPeak), "-", end(masterPeak))
rownames(countMat) = rowRanges_names
head(countMat)

### ==================================== Sequencing depth normalization and differential enriched peaks detection ================ ###
selectR = which(rowSums(countMat) > 5) ## remove low count genes
dataS = countMat[selectR,]
condition = factor(rep(histL, each = length(repL)))
condition
dds = DESeqDataSetFromMatrix(countData = dataS,
                             colData = DataFrame(condition),
                             design = ~ condition)
DDS = DESeq(dds)

normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")

res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")

#Combine results
countMatDiff = cbind(as.data.frame(dataS), as.data.frame(normDDS), as.data.frame(res))
head(countMatDiff)
write.csv(countMatDiff, file = "DESeq2_results_all_counts.csv", row.names = TRUE)
