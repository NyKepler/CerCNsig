#!/usr/bin/env Rscript

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

source("/home/researcher/CerCNsig/Script/main_functions.R")
source("/home/researcher/CerCNsig/Script/helper_functions.R")


num_cores <- 16
type = "HGSC_VS_good"
input <- read.csv(paste0("/home/researcher/CerCNsig/input/", type, ".csv"), header = T, sep = ",", stringsAsFactors = F)

CN_exposure_over_3 <- readRDS("/home/researcher/CerCNsig/Randomization/CN_exposure_over_3.rds")

#' perform filtering on the HGSC_VS_good samples 
#' only 42 HGSC_VS remaining
input <-inner_join(input, CN_exposure_over_3)

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

CN_features <- extractCopynumberFeatures(segment_tables)
VS_components <- readRDS("/home/researcher/CerCNsig/Randomforest_CNsig_filt/HGSC_VS_good_CNsig_over_3.Components.rds")
SxCMat <- generateSampleByComponentMatrix(CN_features, VS_components, cores=1, subcores=num_cores)
pdf(file = "HGSC_VS_good_CNsig_over_3.heatmap.pdf", onefile=FALSE)
NMF::aheatmap(SxCMat, fontsize = 7, Rowv=FALSE, Colv=FALSE, legend = T, breaks=c(seq(0,199,2),500), main="Sample x Component Matrix Heatmap")
dev.off()

# Basis refers to the signature-by-variable matrix, coefficients refers to patient-by-signature matrix, and consensus refers to the connectivity matrix of patients clustered by their dominant signature across 1000 runs. Best-fit is the run that showed the lowest objective score across the 1000 runs. A value of 6 defines the point of stability in the cophenetic, dispersion and silhouette coefficients, and is the maximum sparsity achievable above the null model (randomness) in the Basis plot.

pdf(file = paste0(type, "_CNsig_over_3.numsig.pdf"), onefile=FALSE)
chooseNumberSignatures(SxCMat)
dev.off()


# Same as "generateSignatures" function provides a wrapper to nmf functions required for determining the number of copy-number signatures in a dataset.
Sample_ids<-rownames(SxCMat)
nsig <- 6 #number of signatures can't be smaller than cohort size when calculating coefficiency.
seed <- 77777
nmfalg<-"brunet"
sigs <- NMF::nmf(t(SxCMat), nsig, seed=seed, method=nmfalg, nrun=1000, .opt = paste0("p",num_cores))
saveRDS(sigs, file = paste0(type, "_CNsig_over_3.signatures.rds"))

pdf(file = paste0(type, "_CNsig_over_3.coefmap.pdf"), onefile=FALSE)
coefmap(sigs, Colv="consensus",tracks=c("consensus:"), main="Signature x Sample Matrix Coefmap")#annCol=annCol,
dev.off

pdf(file = paste0(type, "_CNsig_over_3.basismap.pdf"), onefile=FALSE)
basismap(sigs, Rowv=NA, main="Component x Signature Matrix Basismap")
dev.off()

#' Extract component-by-signature matrix
CxSMat <- NMF::basis(sigs) 
colnames(CxSMat) <- paste0("CerCN", 1:7)
saveRDS(CxSMat, "Component_by_Signature_HGSC_VS_good_CNsig_over_3.rds")

#' Extract signature-by-sample matrix
SxSMat <- NMF::scoef(sigs)
rownames(SxSMat)<-paste0("CerCN",1:6)

#Given a sample-by-component matrix this function quantifies signature exposures using the LCD function from the YAPSA package, returning a normalised signature-by-sample matrix. If the component_by_signature matrix is specified then this matrix is used to define the signatures otherwise the signature definitions from the manuscript are used.
#sigs <- readRDS("/home/minerva/WGS/CNsignatures/Data/MaNiLa/CerCNsig_HGSC_all/HGSC.VS.all.signatures.rds")

CerCNsig <- quantifySignatures(SxCMat, CxSMat)
saveRDS(CerCNsig, "CerCNsig.rds")

sessionInfo() %>% capture.output(file="CerCNsig_session_info.txt")

