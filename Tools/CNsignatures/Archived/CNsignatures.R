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

source("/home/minerva/WGS/CNsignatures/Script/main_functions.R")
source("/home/minerva/WGS/CNsignatures/Script/helper_functions.R")

num_cores <- 4
input <- read.csv("/home/minerva/WGS/CNsignatures/Data/MaNiLa/input.final.revised.2.csv", header = T, sep = ",", stringsAsFactors = F)

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



## Extract CN features:
CN_features <- extractCopynumberFeatures(segment_tables)


## Generate sample-by-component matrix using Brenton HGSC compoents
SxCMat <- generateSampleByComponentMatrix(CN_features, cores=1, subcores=num_cores)

## import tuCNsig feature signature matrix 
CxSMat <- readRDS("/home/minerva/WGS/CNsignatures/Data/MaNiLa/feat_sig_mat.rds")
tu_signatures <- quantifySignatures(SxCMat, CxSMat)

## plot signatures histogram
cbPalette <- c(RColorBrewer::brewer.pal(8,"Dark2"),RColorBrewer::brewer.pal(9,"Set1"),"black")
ggplot2::theme_set(ggplot2::theme_gray(base_size = 10))
my_theme<-ggplot2::theme_bw()+ggplot2::theme(axis.text=ggplot2::element_text(size=10),
                   axis.title=ggplot2::element_text(size=10),
                   strip.text.x = ggplot2::element_text(size = 14),
                   strip.text.y = ggplot2::element_text(size = 14),
                   legend.text = ggplot2::element_text(size = 14),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())


signatures <- data.frame(tu_signatures)
input <- data.frame(input)
matched <- input[input$Library %in% colnames(signatures), ]
colnames(signatures) <- matched$Sample
patient <- matched$Patient[1]
signatures$ID <- rownames(signatures)
signatures <- signatures[,order(signatures[nrow(signatures),])]
sig.melt <- melt(signatures, id.var = 'ID')
sig.melt <- within(sig.melt, ID <- factor(ID, c('s1','s2','s3','s4','s5','s6','s7'), ordered = TRUE))

p1<-ggplot(sig.melt, aes(x = variable, y = value, fill = ID)) + 
  geom_bar(stat = 'identity') + xlab("") + ylab("Exposures") + 
  guides(fill=guide_legend(title = "CNsig")) + 
  scale_fill_manual(values=cbPalette) + theme(axis.text.x = element_text(angle = -90))

matched$Sample <-factor(matched$Sample, levels = levels(sig.melt$variable))
table <- data.frame(t(matched[c(3,4,6,7)][order(matched$Sample),])) %>% rownames_to_column()
table[1,] <- gsub("VS", "", table[1,])
table[1,] <- gsub("MaNiLa", "DNAgard", table[1,])
colnames(table) <- 0:(ncol(table)-1)
table<-table[c(1,3,4),] %>% dplyr::rename(ID = 1) %>% dplyr::mutate(ID = recode(ID, Type = 'Sample Origin*', Cellularity = 'Tumor Fraction'))
t <-ggtexttable(table,rows = NULL, theme = ttheme("light"))
text <- "     *'Archival' cervical samples have been preserved by standard clinical precedure and stored at RT or 4 degree while 'DNAgard' cervical sample has been stored in DNAgard buffer at -80 degree."
text.p <- ggparagraph(text = text, face = "italic", size = 11, color = "black")
p2<-ggarrange(t, text.p, ncol = 1, nrow = 2, align = "h", heights = c(0.2, 0.1))

n <- nrow(matched)
if (n>15) {
  pdf(file = "CNsignatures.pdf", onefile = FALSE, width = 15, height = 10)
  print(
    ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1.5, 0.5))
  )
  dev.off() 
} else if (n>7) {
  pdf(file = "CNsignatures.pdf", onefile = FALSE, width = 12, height = 10)
  print(
    ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1.5, 0.5))
  )
  dev.off()
} else {
  pdf(file = "CNsignatures.pdf", onefile = FALSE, width = 1.6*n, height = 10)
  print(
    ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(1.5, 0.5))
  )
  dev.off()
}

save.image(file = "CNsig.Run.rdata") #save the whole R session
