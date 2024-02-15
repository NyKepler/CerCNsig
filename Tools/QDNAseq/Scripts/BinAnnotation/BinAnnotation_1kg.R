#!/usr/bin/env Rscript

#conda activate R

library(base)
library(utils)
library(QDNAseq)
library(Biobase)
library(BSgenome.Hsapiens.UCSC.hg19)
library(R.cache)
library(digest)
future::plan("multicore", workers=4) #Multiprocess with 4 cores
options(future.globals.maxSize = +Inf)


#Create bins, add mappability and blacklist
binSize <- c(1, 5, 10, 15, 30, 50, 100, 500, 1000)
for (i in binSize)
{ bins <- createBins(bsgenome = BSgenome.Hsapiens.UCSC.hg19, binSize = i)
bins$mappability <- calculateMappability(bins, bigWigFile = "wgEncodeCrgMapabilityAlign100mer.bigWig", bigWigAverageOverBed = "/home/minerva/miniconda3/envs/R/bin/bigWigAverageOverBed")
bins$blacklist <- calculateBlacklist(bins, bedFiles = c("wgEncodeDukeMapabilityRegionsExcludable.bed", "wgEncodeDacMapabilityConsensusExcludable.bed"))

#Calculate residuals
#samples <- read.table("samples_aligned.txt", header=F, stringsAsFactors=F)
#bams = paste0("fastq/", samples$V1, ".sort.bam")
ctrl <- binReadCounts(bins, path="BAM/")
ctrl <- applyFilters(ctrl, residual=FALSE, blacklist=FALSE, mappability=FALSE, bases=FALSE)
bins$residual <- iterateResiduals(ctrl)
bins$use <- bins$chromosome %in% as.character(1:22) & bins$bases > 0

#Convert to AnnotatedDataFrame and add metadata
bins <- AnnotatedDataFrame(bins, varMetadata = data.frame(labelDescription=c(
"Chromosome name",
"Base pair start position",
"Base pair end position",
"Percentage of non-N nucleotides (of full bin size)",
"Percentage of C and G nucleotides (of non-N nucleotides)",
"Average mappability of 100mers with a maximum of 2 mismatches",
"Percent overlap with ENCODE blacklisted regions",
"Median loess residual from 1000 Genomes PE LPWGS (2x100)",
"Whether the bin should be used in subsequent analysis steps"), row.names = colnames(bins)))

#Additional descriptive metadata 
attr(bins, paste0("h19.PE100.", i, "kb")) <- list(
author="Minerva Li",
date=Sys.time(),
organism="Hsapiens",
build="hg19",
version="0.0.2",
md5=digest::digest(bins@data),
sessionInfo=sessionInfo())


#Compare to pre-calculatedAnnotations
bins_precalc <- getBinAnnotations(binSize=i)
sum((pData(bins_precalc)$bases - pData(bins)$bases))
sum((pData(bins_precalc)$gc - pData(bins)$gc), na.rm=TRUE)
sum((pData(bins_precalc)$mappability - pData(bins)$mappability))
sum((pData(bins_precalc)$blacklist - pData(bins)$blacklist))
cor.test(pData(bins_precalc)$residual, pData(bins)$residual)
sum(pData(bins_precalc)$use)
sum(pData(bins)$use)

saveRDS(bins, file = paste0("BinAnnotation/hg19_PE100.", i, "kb.rds"))}




