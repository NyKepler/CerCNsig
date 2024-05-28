## Load library
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
library(RColorBrewer)
library(scales)

my_theme<-ggplot2::theme_bw()+ggplot2::theme(axis.text=ggplot2::element_text(size=10),
                                             axis.title=ggplot2::element_text(size=10),
                                             strip.text.x = ggplot2::element_text(size = 14),
                                             strip.text.y = ggplot2::element_text(size = 14),
                                             legend.text = ggplot2::element_text(size = 14),
                                             panel.grid.minor = ggplot2::element_blank(),
                                             panel.grid.major = ggplot2::element_blank())

## 1. Import data
all.sample <- read.csv("/home/researcher/CerCNsig/input/All_Samples.csv", header = T, sep = ",", stringsAsFactors = F)
folder.name <- "Randomforest_CNsig_filt/CerCNsig_6"
setwd(paste(path_home_r(), "CerCNsig", folder.name, sep = "/"))
SSMatrix.t <- readRDS("CerCNsig.Tissue_Tumor.SSmatrix.rds")
SSMatrix.vs <- readRDS("CerCNsig.All_Cervical.SSmatrix.rds")

## Tissue and tumor samples
#' Subset different group
SSMatrix <- SSMatrix.t
SSMatrix <- left_join(SSMatrix, all.sample, by = join_by(sample == Library))
SSMatrix.HGSC.t <- filter(SSMatrix, SSMatrix$Group == "HGSC" & SSMatrix$Type == "Tumor")
SSMatrix.HGSC.ti <- filter(SSMatrix, SSMatrix$Group == "HGSC" & SSMatrix$Type == "Normal tissues") %>% mutate(Group = "HGSC_NT")
SSMatrix.RRSO <- filter(SSMatrix, SSMatrix$Group == "BRCA_RRSO") %>% mutate(Group = "RRSO")
SSMatrix.Benign <- filter(SSMatrix, SSMatrix$Group == "Benign")


SSMatrix.HGSC.t <-column_to_rownames(SSMatrix.HGSC.t, var = "sample") %>% select(paste0("CerCN", 1:6), "Type", "Group") %>% reshape2::melt() 
SSMatrix.HGSC.ti <-column_to_rownames(SSMatrix.HGSC.ti, var = "sample") %>% select(paste0("CerCN", 1:6), "Type", "Group") %>% reshape2::melt() 
SSMatrix.RRSO <-column_to_rownames(SSMatrix.RRSO, var = "sample") %>% select(paste0("CerCN", 1:6), "Type", "Group") %>% reshape2::melt() 
SSMatrix.Benign <-column_to_rownames(SSMatrix.Benign, var = "sample") %>% select(paste0("CerCN", 1:6), "Type", "Group") %>% reshape2::melt() 

#' skip the normal tissue from HGSC 
pdat_tf<-rbind(SSMatrix.Benign, SSMatrix.RRSO, SSMatrix.HGSC.t)
colnames(pdat_tf) <- c("Type", "Group", "Signature", "Cosine")
pdat_tf$Group <- factor(pdat_tf$Group, levels = c("Benign", "RRSO", "HGSC"))
signif_tf <-ggpubr::compare_means(Cosine ~ Group,
                                  data=pdat_tf,
                                  group.by = c("Signature"),
                                  method="kruskal.test",
                                  p.adjust.method="BH")

signif_tf<-signif_tf %>% filter(p.signif!="ns") %>% group_by(Signature) %>% filter(p.format==min(p.format))
signif_tf$Cosine<-1
knitr::kable(signif_tf)
write.table(signif_tf, "CerCNsig_Cosine_Comparison_Tissue_Tumor.csv", row.names = F,
            col.names = T, quote = F, sep = ",")

p.tf <- ggplot(pdat_tf,aes(x=Signature,y=Cosine,fill=Group))+
  geom_boxplot()+
  my_theme+
  scale_fill_manual(values = c("#8DD3C7", "#FFFFB3", "#FB8072")) + 
  coord_cartesian(ylim=c(0,1))+
  annotate("text",x=signif_tf$Signature,y=signif_tf$Cosine,label=signif_tf$p.signif, size = 6) + 
  theme(legend.position="right",legend.key.size = unit(0.4, 'cm'), legend.title = element_blank())
ggsave("CerCNsig_Cosine_Comparison_Tissue_Tumor.pdf", width = 6.13, height = 5.16, units = "in", dpi = 600)

#' include the normal tissue from HGSC
pdat_tf<-rbind(SSMatrix.Benign, SSMatrix.RRSO, SSMatrix.HGSC.t, SSMatrix.HGSC.ti)
colnames(pdat_tf) <- c("Type", "Group", "Signature", "Cosine")
pdat_tf$Group <- factor(pdat_tf$Group, levels = c("Benign", "RRSO", "HGSC_NT", "HGSC"))
signif_tf <-ggpubr::compare_means(Cosine ~ Group,
                                  data=pdat_tf,
                                  group.by = c("Signature"),
                                  method="kruskal.test",
                                  p.adjust.method="BH")

signif_tf<-signif_tf %>% filter(p.signif!="ns") %>% group_by(Signature) %>% filter(p.format==min(p.format))
signif_tf$Cosine<-1
knitr::kable(signif_tf)
write.table(signif_tf, "CerCNsig_Cosine_Comparison_Tissue_Tumor_2.csv", row.names = F,
            col.names = T, quote = F, sep = ",")

p.tf <- ggplot(pdat_tf,aes(x=Signature,y=Cosine,fill=Group))+
  geom_boxplot()+
  my_theme+
  scale_fill_manual(values = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072")) + 
  coord_cartesian(ylim=c(0,1))+
  annotate("text",x=signif_tf$Signature,y=signif_tf$Cosine,label=signif_tf$p.signif, size = 6) + 
  theme(legend.position="right",legend.key.size = unit(0.4, 'cm'), legend.title = element_blank())
ggsave("CerCNsig_Cosine_Comparison_Tissue_Tumor_2.pdf", width = 6.13, height = 5.16, units = "in", dpi = 600)

#' Only normal tissue vs tumor in HGSC
pdat_tf<-rbind(SSMatrix.HGSC.t, SSMatrix.HGSC.ti)
colnames(pdat_tf) <- c("Type", "Group", "Signature", "Cosine")
pdat_tf$Group <- factor(pdat_tf$Group, levels = c("HGSC_NT", "HGSC"))
signif_tf <-ggpubr::compare_means(Cosine ~ Group,
                                  data=pdat_tf,
                                  group.by = c("Signature"),
                                  method="wilcox.test",
                                  p.adjust.method="BH")

signif_tf<-signif_tf %>% filter(p.signif!="ns") %>% group_by(Signature) %>% filter(p.format==min(p.format))
signif_tf$Cosine<-1
knitr::kable(signif_tf)
write.table(signif_tf, "CerCNsig_Cosine_Comparison_Tissue_Tumor_HGSC.csv", row.names = F,
            col.names = T, quote = F, sep = ",")

p.tf <- ggplot(pdat_tf,aes(x=Signature,y=Cosine,fill=Group))+
  geom_boxplot()+
  my_theme+
  scale_fill_manual(values = c("#BEBADA", "#FB8072")) + 
  coord_cartesian(ylim=c(0,1))+
  annotate("text",x=signif_tf$Signature,y=signif_tf$Cosine,label=signif_tf$p.signif, size = 6) + 
  theme(legend.position="right",legend.key.size = unit(0.4, 'cm'), legend.title = element_blank())
ggsave("CerCNsig_Cosine_Comparison_Tissue_Tumor_HGSC.pdf", width = 6.13, height = 5.16, units = "in", dpi = 600)

## Cervical samples
#' Empty enviroment before running this
#' Subset different group
SSMatrix <- SSMatrix.vs
SSMatrix <- left_join(SSMatrix, all.sample, by = join_by(sample == Library))
SSMatrix.HGSC <- filter(SSMatrix, SSMatrix$Group == "HGSC")
SSMatrix.RRSO <- filter(SSMatrix, SSMatrix$Group == "BRCA_RRSO") %>% mutate(Group = "RRSO")
SSMatrix.Benign <- filter(SSMatrix, SSMatrix$Group == "Benign")


SSMatrix.HGSC <-column_to_rownames(SSMatrix.HGSC, var = "sample") %>% select(paste0("CerCN", 1:6), "Type", "Group") %>% reshape2::melt() 
SSMatrix.RRSO <-column_to_rownames(SSMatrix.RRSO, var = "sample") %>% select(paste0("CerCN", 1:6), "Type", "Group") %>% reshape2::melt() 
SSMatrix.Benign <-column_to_rownames(SSMatrix.Benign, var = "sample") %>% select(paste0("CerCN", 1:6), "Type", "Group") %>% reshape2::melt() 

#' compare VS in 3 groups
pdat_tf<-rbind(SSMatrix.Benign, SSMatrix.RRSO, SSMatrix.HGSC)
colnames(pdat_tf) <- c("Type", "Group", "Signature", "Cosine")
pdat_tf$Group <- factor(pdat_tf$Group, levels = c("Benign", "RRSO", "HGSC"))
signif_tf <-ggpubr::compare_means(Cosine ~ Group,
                                  data=pdat_tf,
                                  group.by = c("Signature"),
                                  method="kruskal.test",
                                  p.adjust.method="BH")

signif_tf<-signif_tf %>% filter(p.signif!="ns") %>% group_by(Signature) %>% filter(p.format==min(p.format))
signif_tf$Cosine<-1
knitr::kable(signif_tf)
write.table(signif_tf, "CerCNsig_Cosine_Comparison_VS.csv", row.names = F,
            col.names = T, quote = F, sep = ",")

p.tf <- ggplot(pdat_tf,aes(x=Signature,y=Cosine,fill=Group))+
  geom_boxplot()+
  my_theme+
  scale_fill_manual(values = c("#8DD3C7", "#FFFFB3", "#FB8072")) + 
  coord_cartesian(ylim=c(0,1))+
  annotate("text",x=signif_tf$Signature,y=signif_tf$Cosine,label=signif_tf$p.signif, size = 6) + 
  theme(legend.position="right",legend.key.size = unit(0.4, 'cm'), legend.title = element_blank())
ggsave("CerCNsig_Cosine_Comparison_VS.pdf", width = 6.13, height = 5.16, units = "in", dpi = 600)

#' compare only the VS in HGSC group
SSMatrix.HGSC <- filter(SSMatrix, SSMatrix$Group == "HGSC")
SSMatrix.HGSC.pdo6 <- filter(SSMatrix.HGSC, SSMatrix.HGSC$Sample.group == "HGSC_VS_prediagnostic_>6") %>% mutate(Sample.group = "Prediag_>6m")
SSMatrix.HGSC.pdu6 <- filter(SSMatrix.HGSC, SSMatrix.HGSC$Sample.group == "HGSC_VS_prediagnostic_0-6") %>% mutate(Sample.group = "Prediag_0-6m")
SSMatrix.HGSC.d <- filter(SSMatrix.HGSC, SSMatrix.HGSC$Sample.group == "HGSC_VS_diagnostic") %>% mutate(Sample.group = "Diagnostic")

SSMatrix.HGSC.pdo6 <-column_to_rownames(SSMatrix.HGSC.pdo6, var = "sample") %>% select(paste0("CerCN", 1:6), "Type", "Sample.group") %>% reshape2::melt() 
SSMatrix.HGSC.pdu6 <-column_to_rownames(SSMatrix.HGSC.pdu6, var = "sample") %>% select(paste0("CerCN", 1:6), "Type", "Sample.group") %>% reshape2::melt() 
SSMatrix.HGSC.d <-column_to_rownames(SSMatrix.HGSC.d, var = "sample") %>% select(paste0("CerCN", 1:6), "Type", "Sample.group") %>% reshape2::melt() 

pdat_tf<-rbind(SSMatrix.HGSC.pdo6, SSMatrix.HGSC.pdu6, SSMatrix.HGSC.d)
colnames(pdat_tf) <- c("Type", "Group", "Signature", "Cosine")
pdat_tf$Group <- factor(pdat_tf$Group, levels = c("Prediag_>6m", 
                                                  "Prediag_0-6m", 
                                                  "Diagnostic"))
signif_tf <-ggpubr::compare_means(Cosine ~ Group,
                                  data=pdat_tf,
                                  group.by = c("Signature"),
                                  method="kruskal.test",
                                  p.adjust.method="BH")

signif_tf<-signif_tf %>% filter(p.signif!="ns") %>% group_by(Signature) %>% filter(p.format==min(p.format))
signif_tf$Cosine<-1
knitr::kable(signif_tf)
write.table(signif_tf, "CerCNsig_Cosine_Comparison_VS_HGSC.csv", row.names = F,
            col.names = T, quote = F, sep = ",")

p.tf <- ggplot(pdat_tf,aes(x=Signature,y=Cosine,fill=Group))+
  geom_boxplot()+
  my_theme+
  scale_fill_manual(values = c("#8DD3C7", "#FFFFB3", "#FB8072")) + 
  coord_cartesian(ylim=c(0,1))+
  annotate("text",x=signif_tf$Signature,y=signif_tf$Cosine,label=signif_tf$p.signif, size = 6) + 
  theme(legend.position="right",legend.key.size = unit(0.4, 'cm'), legend.title = element_blank())
ggsave("CerCNsig_Cosine_Comparison_VS_HGSC.pdf", width = 6.8, height = 5.16, units = "in", dpi = 600)