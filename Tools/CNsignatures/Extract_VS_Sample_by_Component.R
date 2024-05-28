#!/usr/bin/env Rscript

library(dplyr)
library(flexmix)
library(QDNAseq)
library(YAPSA)
library(stringr)
library(ggplot2)
library(ggpubr)
library(reshape)
library(tidyverse)
library(tibble)
library(patchwork)

source("/home/researcher/CerCNsig/Script/main_functions.R")
source("/home/researcher/CerCNsig/Script/helper_functions.R")

num_cores <- 16
type = "Normal_Tissue"
input <- read.csv(paste0("/home/researcher/CerCNsig/input/", type, ".csv"), header = T, sep = ",", stringsAsFactors = F)

set.seed(1234)
# collapse equal value neighboring segments
segment_tables <-list()
for(i in input$Library)
{
  if(file.exists(paste0(i, "_seg.tsv")))
  {
    cn <- read.table(paste0(i, "_seg.tsv"),header=T,
                   sep="\t",stringsAsFactors = F)             	
    if (nrow(cn)>0)
    {	
    segment_tables[[i]] <- cn
    }
  }
}

saveRDS(segment_tables, file = paste0(type, ".Segtable.rds"))


# This function takes as input a collection of absolute copy-number profiles and returns a list of copy-number features extracted from these samples. Six fundamental CN features: the number of breakpoints per 10MB, the copy-number of segments, the size of segments, the difference in CN between adjacent segments, the distance of breakpoints from the centromere, and the lengths of oscillating CN segment chains.
CN_features <- extractCopynumberFeatures(segment_tables)
saveRDS(CN_features, file = paste0(type, ".CNfeatures.rds"))

# Optional! This function takes as input a list of copy number features and outputs a list of mixture model components fit to these features.
# VS_components <- fitMixtureModels(CN_features)
# saveRDS(VS_components, file = paste0(type, ".Components.rds"))

VS_components <- readRDS("/home/researcher/CerCNsig/Randomforest_CNsig/HGSC_VS_good.Components.rds")

# Given a set of extracted copy number features profiles this function returns a sum-of-posterior sample-by-component matrix. If the all_components argument is specified, then the sum-of-posteriors is calculated using these components, otherwise the component definitions from the manuscript are used.
sample_component_matrix <- generateSampleByComponentMatrix(CN_features, VS_components, cores=1, subcores=num_cores)
pdf(file = paste0(type, ".heatmap.pdf"), onefile=FALSE)
NMF::aheatmap(sample_component_matrix, fontsize = 7, Rowv=FALSE, Colv=FALSE, legend = T, breaks=c(seq(0,199,2),500), main="Component x Sample matrix Heatmap")
dev.off()
saveRDS(sample_component_matrix, file = paste0(type, ".SxCMat.rds"))

save.image(file = paste0(type, ".rdata")) #save the whole R session

sessionInfo() %>% capture.output(file="session_info.txt")
