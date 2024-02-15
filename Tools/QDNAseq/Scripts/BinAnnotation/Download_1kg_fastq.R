#!/usr/bin/env Rscript

#conda activate R

library(base)
library(utils)
library(stats)
library(plyr)


# Downloading 1000 Genomes samples
###################################################
## # download table of samples
urlroot <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3"
g1k <- read.table(file.path(urlroot, "20130502.phase3.sequence.index"), header=TRUE, sep="\t", as.is=TRUE, fill=TRUE)

# Keep cases that are on Illumina platform, low coverage, not withdrawn and paired-end sequenced in a single lane.
g1k <- g1k[g1k$INSTRUMENT_PLATFORM == "ILLUMINA", ]
g1k <- g1k[g1k$INSTRUMENT_MODEL == "Illumina HiSeq 2000", ]
g1k <- g1k[g1k$ANALYSIS_GROUP == "low coverage", ]
g1k <- g1k[g1k$LIBRARY_LAYOUT == "PAIRED", ]
g1k <- g1k[g1k$WITHDRAWN == 0, ]
g1k <- g1k[!g1k$PAIRED_FASTQ == "", ]

# Keep cases with total read lengths of at least 200bp, library insert size not bigger than 500bp and total read count below 100M. 
g1k$BASE_COUNT <- as.numeric(g1k$BASE_COUNT)
g1k$READ_COUNT <- as.integer(g1k$READ_COUNT)
g1k$READ_LENGTH <- g1k$BASE_COUNT / g1k$READ_COUNT
readLengthPerRun <- aggregate(g1k$READ_LENGTH, by=list(sample=g1k$RUN_ID), FUN=sum)
g1k <- g1k[g1k$RUN_ID %in% readLengthPerRun$sample[readLengthPerRun$x == 200], ]
g1k <- g1k[g1k$INSERT_SIZE > 0 & g1k$INSERT_SIZE <= 350, ]
readCountPerRun <- aggregate(g1k$READ_COUNT, by=list(sample=g1k$RUN_ID), FUN=sum)
g1k <- g1k[g1k$RUN_ID %in% readCountPerRun$sample[readCountPerRun$x <= 10e+07], ]


# Sample 100 sequencing runs and then 50 individuals which resulting 38 paired-end sequenced genomes from 15 different population.
set.seed(7777)
g1k <- g1k[g1k$POPULATION %in% sample(g1k$POPULATION, 100), ]
g1k <- g1k[g1k$RUN_ID %in% sample(g1k$RUN_ID, 100), ]
freq <- count(g1k, 'SAMPLE_ID')
g1k<- g1k[g1k$SAMPLE_ID %in% freq$SAMPLE_ID[freq$freq == 2], ]
g1k$fileName <- basename(g1k$FASTQ_FILE)

# Download FASTQ files
seqroot = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3" #latest data 
for (i in rownames(g1k)) { 
  sourceFile <- file.path(seqroot, g1k[i, "FASTQ_FILE"])
  destFile <- g1k[i, "fileName"]
  if (!file.exists(destFile))
    download.file(sourceFile, destFile, method = "wget")
    
}

x = unique(subset(g1k, select=c(SAMPLE_NAME, RUN_ID)))
write.table(x, "sample_fastq_map.csv", col.names=T, quote=F, row.names=F, sep=",")

