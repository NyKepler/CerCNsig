#!/usr/bin/env Rscript

library(base)
library(utils)
library(R.cache)
library(digest)
library(stringr)
library(QDNAseq)
library(DNAcopy)
library(CGHcall)
library(CGHregions)
library(Biobase)
future::plan("multisession", workers=4) #Multiprocess with 4 cores
options(future.globals.maxSize = +Inf)
options("QDNAseq::verbose" = NA) #Send verbose output to stdout instead of stderr


# Extract sample info and bin annotation file
input <- read.csv("/home/minerva/WGS/QDNAseq/Data/MaNiLa/MaNiLa_CS3_input.csv", header=T, sep=",", stringsAsFactors = F)
sampleID = str_sub(dir(pattern = ".bam.bai"), end = -9)
binsize <- input[input$library %in% sampleID, ]$bin

# QDNAseq bin annotation for hg19: pre-calculated for genome buildhg19 and bin sizes 1, 5, 10, 15, 30, 50, 100, 500, and 1000 kbp, otherwise use createBins(). For deeper sequencing depths (up to 1x) smaller binsize 5kb gives greater resolution, 0.1-0.15X depthy 15kb is optimal and up to 1X 5kb is ok.
bins <- readRDS(file=paste0("/home/minerva/WGS/QDNAseq/BinAnnotation_v2/BinAnnotation/hg19_PE100.", binsize, "kb.rds"))

# Processing BAM file to extract read counts
readCounts <- binReadCounts(bins, cache = TRUE, chunkSize = 1E7, pairedEnds = TRUE)  
main = sampleID
png(file = paste0(sampleID, ".", binsize, "kb.rc.png"), width = 1000, height = 500)
plot(readCounts, main = main, logTransform = FALSE, ylim = c(-50,200))   # or continue with ReadcountPlot R script to plot with for loops
highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE) 
dev.off()

# Export read counts for ichorCNA or other analysis
saveRDS(readCounts, file = paste0(sampleID, ".", binsize, "kb.rc.rds")) 


# Filter and normalize reads
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
png(file = paste0(sampleID, ".", binsize, "kb.isobar.png"), width = 1000, height = 500)
isobarPlot(readCountsFiltered, main = main)
dev.off()

readCountsFiltered <- estimateCorrection(readCountsFiltered)
png(file = paste0(sampleID, ".", binsize, "kb.noise.png"), width = 1000, height = 500)
noisePlot(readCountsFiltered, main = main)
dev.off()

copyNumbers <- correctBins(readCountsFiltered) #Correction for GC content and mappabilit
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized) #Data will be log2 transfer in this step. 
png(file = paste0(sampleID, ".", binsize, "kb.sm.png"), width = 1000, height = 500)
plot(copyNumbersSmooth, main = main)
dev.off()


# Downstream analysis with DNAcopy and CGHcall: 
# 1. Copy Number Segmentation
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun = "sqrt") #log2-transformation or a sqrt(x + 3/8) which stabilizes the variance of a Poisson distribution (Anscombe transform). For downstream analysis using ABSOLUTE, choose "none" in the transformFun.
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
png(file = paste0(sampleID, ".", binsize, "kb.sg.png"), width = 1000, height = 500)
plot(copyNumbersSegmented, main = main)
dev.off()
saveRDS(copyNumbersSegmented, file = paste0(sampleID, ".", binsize, "kb.sg.rds")) 

# 2. Copy Number calling
copyNumbersCalled <- callBins(copyNumbersSegmented, organism="human", method="CGHcall")
png(file = paste0(sampleID, ".", binsize, "kb.cn.png"), width = 1000, height = 500)
plot(copyNumbersCalled, main = main)
dev.off()

png(file = paste0(sampleID, ".", binsize, "kb.freq.png"), width = 1000, height = 500)
frequencyPlot(copyNumbersCalled, main = main) # or use Copynumberscall R script to plot for all samples and for all chromosomes seperately, otherwise
dev.off()


# In our experience CNA are most clearly depicted by loading both the Copynumber.igv and the Segmentation.igv files concurrently, then changing the type of graph to a bar chart for the segmentation track by right clicking the track and selecting bar chart.
exportBins(copyNumbersCalled, file=paste0(sampleID, ".", binsize, "kb.sg.igv"), format="igv", type="segments")
exportBins(copyNumbersCalled, file=paste0(sampleID, ".", binsize, "kb.cn.igv"), format="igv", type="copynumber")

saveRDS(copyNumbersCalled, file = paste0(sampleID, ".", binsize, "kb.cn.rds")) 


cgh <- makeCgh(copyNumbersCalled) #if using CGHregions for downstream analysis


