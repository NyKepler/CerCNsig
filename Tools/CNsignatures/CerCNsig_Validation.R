#' Calculate cosine based on Guyuan 
#' BINP52_CNA_Framework/workflow/scripts/CN_sig.R
library(NMF)
library(flexmix)
library(YAPSA)
library(tidyverse)
library(lsa)

# apply cosine similarity to find the enriched signature for each sample
## 1. Import data
folder.name <- "Randomforest_CNsig_filt"
setwd(paste(path_home_r(), "CerCNsig", folder.name, sep = "/"))

#' read RDS files and combine benign and RRSO group into one group called 'Benign'
Benign_VS_good.SxCMat <- readRDS("Benign_VS_good.SxCMat.rds") %>% as.data.frame()
Benign_VS_good.SxCMat$Class <- "Benign"
HGSC_VS_good.SxCMat <- readRDS("HGSC_VS_good.SxCMat.rds") %>% as.data.frame()
HGSC_VS_good.SxCMat$Class <- "HGSC"
RRSO_VS.SxCMat <- readRDS("RRSO_VS.SxCMat.rds") %>% as.data.frame()
RRSO_VS.SxCMat$Class <- "RRSO"

VS <-rbind(Benign_VS_good.SxCMat, HGSC_VS_good.SxCMat, RRSO_VS.SxCMat)

CxSMat <- readRDS("Component_by_Signature_HGSC_VS_good_CNsig_over_3.rds")
CerCNsig <- paste0('CerCN',c(1:7))

## load the Sample-by-Component matrix generated above
SCmatrix <- t(VS[1:36])
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
SSmatrix <- SSmatrix %>% select('sample', 'enrich', 'CerCN1':'CerCN7')
row.names(SSmatrix) <- 1:nrow(SSmatrix)
write.table(SSmatrix, file = "CerCNsig.good.VS.SSmatrix.tsv", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
saveRDS(SSmatrix, file = "CerCNsig.good.VS.SSmatrix.rds", compress = FALSE)