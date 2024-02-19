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

source("/home/minerva/WGS/CNsignatures/Script/main_functions.R")
source("/home/minerva/WGS/CNsignatures/Script/helper_functions.R")

input <- read.csv("/home/minerva/WGS/CNsignatures/Data/MaNiLa/input.final.revised.csv", header = T, sep = ",", stringsAsFactors = F)

num_cores <- 4

type = "BRCA.27.CNsig"

set.seed(1234)

# collapse equal value neighboring segments
segment_tables <-list()
for(i in input$Library)
{
  if(file.exists(paste0(i, "_segments.tsv")))
  {
    cn <- read.table(paste0(i, "_segments.tsv"),header=T,
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
#MaNiLa_components <- fitMixtureModels(CN_features)
#saveRDS(MaNiLa_components, file = "/home/minerva/WGS/CNsignatures/Data/MaNiLa/MaNiLa.Tumor.Components.rds")


# Given a set of extracted copy number features profiles this function returns a sum-of-posterior sample-by-component matrix. If the all_components argument is specified, then the sum-of-posteriors is calculated using these components, otherwise the component definitions from the manuscript are used.
sample_component_matrix <- generateSampleByComponentMatrix(CN_features, cores=1, subcores=num_cores)
pdf(file = paste0(type, ".heatmap.pdf"), onefile=FALSE)
NMF::aheatmap(sample_component_matrix, fontsize = 7, Rowv=FALSE, Colv=FALSE, legend = T, breaks=c(seq(0,199,2),500), main="Component x Sample matrix Heatmap")
dev.off()
saveRDS(sample_component_matrix, file = paste0(type, ".SxCMat.rds"))


# Basis refers to the signature-by-variable matrix, coefficients refers to patient-by-signature matrix, and consensus refers to the connectivity matrix of patients clustered by their dominant signature across 1000 runs. Best-fit is the run that showed the lowest objective score across the 1000 runs. A value of 6 defines the point of stability in the cophenetic, dispersion and silhouette coefficients, and is the maximum sparsity achievable above the null model (randomness) in the Basis plot.

pdf(file = paste0(type, ".numsig.pdf"), onefile=FALSE)
chooseNumberSignatures(sample_component_matrix)
dev.off()


# Same as "generateSignatures" function provides a wrapper to nmf functions required for determining the number of copy-number signatures in a dataset.
Sample_ids<-rownames(sample_component_matrix)
nsig <- 7 #number of signatures can't be smaller than cohort size when calculating coefficiency.
seed <- 77777
nmfalg<-"brunet"
sigs <- NMF::nmf(t(sample_component_matrix), nsig, seed=seed, method=nmfalg, nrun=1000, .opt = paste0("p",num_cores))
saveRDS(sigs, file = paste0(type, ".signatures.rds"))

pdf(file = paste0(type, ".coefmap.pdf"), onefile=FALSE)
coefmap(sigs, Colv="consensus",tracks=c("consensus:"), main="Sample x Signature matrix Coefmap")
dev.off

pdf(file = paste0(type, ".basismap.pdf"), onefile=FALSE)
basismap(sigs, Rowv=NA, main="Signature x Component matrix Basismap")
dev.off()

#Given a sample-by-component matrix this function quantifies signature exposures using the LCD function from the YAPSA package, returning a normalised signature-by-sample matrix. If the component_by_signature matrix is specified then this matrix is used to define the signatures otherwise the signature definitions from the manuscript are used.
## apply the MaNiLa tumor signatures instead

sig_thresh<-0.01
sig_pat_mat <- scoef(sigs)
rownames(sig_pat_mat)<-paste0("s",1:nsig)
CN_signatures <-normaliseMatrix(sig_pat_mat,sig_thresh)

cbPalette <- c(RColorBrewer::brewer.pal(8,"Dark2"),RColorBrewer::brewer.pal(9,"Set3"),"black")
ggplot2::theme_set(ggplot2::theme_gray(base_size = 10))
my_theme<-ggplot2::theme_bw()+ggplot2::theme(axis.text=ggplot2::element_text(size=10),
                   axis.title=ggplot2::element_text(size=10),
                   strip.text.x = ggplot2::element_text(size = 14),
                   strip.text.y = ggplot2::element_text(size = 14),
                   legend.text = ggplot2::element_text(size = 14),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())
signatures <- data.frame(CN_signatures)
input <- data.frame(input)
matched <- input[input$Library %in% colnames(signatures), ]
colnames(signatures) <- matched$Sample
signatures$ID <- rownames(signatures)
signatures <- signatures[,order(signatures[nrow(signatures),])]
sig.melt <- melt(signatures, id.var = 'ID')
sig.melt <- within(sig.melt, ID <- factor(ID, c('s1','s2','s3','s4','s5','s6','s7'), ordered = TRUE))

p<-ggplot(sig.melt, aes(x = variable, y = value, fill = ID)) + 
  geom_bar(stat = 'identity') + xlab("") + ylab("Exposures") + ylim(0, 1) +
  guides(fill=guide_legend(title = "CNsig")) + scale_fill_manual(values=cbPalette) + 
  theme(axis.text.x = element_text(angle = -90)) 

matched$Sample <-factor(matched$Sample, levels = levels(sig.melt$variable))
table <- data.frame(t(matched[c(3,4,6,7)][order(matched$Sample),])) %>% rownames_to_column()
colnames(table) <- gsub("HGSOC_", "", table[2,])
table<-table[c(1,3,4),] %>% dplyr::rename(ID = Sample) %>% mutate(ID = recode(ID, Type = 'Sample Origin', Cellularity = 'Tumor Fraction'))
t <-ggtexttable(table,rows = NULL, theme = ttheme("light"))

n <- nrow(matched)
pdf(file = paste0(type, ".Signature.pdf"), onefile = FALSE, width = 12, height = 10)
ggarrange(p, t, ncol = 1, nrow = 2, heights = c(1, 0.5))
dev.off() 

saveRDS(CN_signatures, file = paste0(type, ".CNsignatures.rds")) 
save.image(file = paste0(type, ".rdata")) #save the whole R session
