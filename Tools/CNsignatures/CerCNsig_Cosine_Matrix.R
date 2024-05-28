#' Calculate cosine based on Guyuan 
#' BINP52_CNA_Framework/workflow/scripts/CN_sig.R

library(fs)
library(dplyr)
library(NMF)
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
library(cowplot)
library(lsa)

# apply cosine similarity to find the enriched signature for each sample
## 1. Import data
folder.name <- "Randomforest_CNsig_filt/CerCNsig_6"
setwd(paste(path_home_r(), "CerCNsig", folder.name, sep = "/"))

#' read RDS files and combine benign and RRSO group into one group called 'Benign'
Benign_VS_good.SxCMat <- readRDS("Benign_VS_good.SxCMat.rds") %>% as.data.frame()
Benign_VS_good.SxCMat$Class <- "Benign_VS"
RRSO_VS.SxCMat <- readRDS("RRSO_VS.SxCMat.rds") %>% as.data.frame()
RRSO_VS.SxCMat$Class <- "RRSO_VS"
Normal_Tissue.SxCMat <- readRDS("Normal_Tissue.SxCMat.rds") %>% as.data.frame()
Normal_Tissue.SxCMat$Class <- "Normal_Tissue"
HGSC_VS_good.SxCMat <- readRDS("HGSC_VS_good.SxCMat.rds") %>% as.data.frame()
HGSC_VS_good.SxCMat$Class <- "HGSC_VS"
HGSC_Tumor.SxCMat <- readRDS("HGSC_Tumor.SxCMat.rds") %>% as.data.frame()
HGSC_Tumor.SxCMat$Class <- "HGSC_Tumor"
SxCMat <-rbind(Normal_Tissue.SxCMat, HGSC_Tumor.SxCMat)

#' All Cervical samples
Cervical.SxCMat <- readRDS("All_Cervical.SxCMat.rds") %>% as.data.frame()
Cervical.SxCMat$Class <- "Cervical"
SxCMat <- Cervical.SxCMat
name <- "Tissue_Tumor" #' Tissue_Tumor_VS, all_cervical


CxSMat <- readRDS("Component_by_Signature_HGSC_VS_good_CNsig_over_3.rds")
CerCNsig <- paste0('CerCN',c(1:6))

## load the Sample-by-Component matrix generated above
SCmatrix <- t(SxCMat[1:36])
samples <- colnames(SCmatrix)

## calculate the cosine similarity
cos_mat <- cbind(CxSMat, SCmatrix)
cos_mat <- cosine(cos_mat)

## extract the required signature similarity matrix
SSmatrix <- cos_mat[samples,]
SSmatrix <- as.data.frame(SSmatrix) %>% select(all_of(CerCNsig))

## find the most similar signature for each sample
sample_by_signature <- apply(SSmatrix, 1, function(t) colnames(SSmatrix)[which.max(t)]) 
sample_by_signature <- as.list(sample_by_signature)
SSmatrix['enrich'] <- NA
for (i in samples) {
  sigID <- sample_by_signature[[i]]
  SSmatrix[i, 'enrich'] <- sigID
}

## save the SSmatrix
SSmatrix['sample'] <- rownames(SSmatrix)
SSmatrix <- SSmatrix %>% select('sample', 'enrich', 'CerCN1':'CerCN6')
row.names(SSmatrix) <- 1:nrow(SSmatrix)
write.table(SSmatrix, file = paste0("CerCNsig.", name, ".SSmatrix.tsv"), sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
saveRDS(SSmatrix, file = paste0("CerCNsig.", name, ".SSmatrix.rds"), compress = FALSE)

sessionInfo() %>% capture.output(file="CerCNsig_Cosine_session_info.txt")
