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

## 1. Import data
folder.name <- "Randomforest_CNsig_filt/CerCNsig_6"
setwd(paste(path_home_r(), "CerCNsig", folder.name, sep = "/"))

#' import component-by-signature matrix from CerCNsig
CxSMat <- readRDS("Component_by_Signature_HGSC_VS_good_CNsig_over_3.rds")
SxCMat <- t(CxSMat)

temp<-as.data.frame(SxCMat) 
norm_const<-apply(temp,2,sum)
temp<-data.frame(t(apply(temp,1,function(x){x/norm_const})))
temp$sig<-rownames(SxCMat)
pdat<-reshape2::melt(temp,id.vars="sig")
vars<-gsub("^\\d+|\\d+$", "",as.character(pdat$variable))
pdat<-cbind(pdat,vars)
colnames(pdat)<-c("sig","Feature","value","Distribution")
pdat$sig<-factor(pdat$sig,levels=paste0("CerCN", 1:6))

pdat$Distribution<-plyr::revalue(pdat$Distribution,
                                 c(bp10MB="Breakpoint number",
                                   copynumber="Copy-number",
                                   changepoint="CN changepoint",
                                   bpchrarm="Breakpoints per chr arm",
                                   osCN="Oscilating CN length",
                                   segsize="Segment size"))
pdat$Distribution<-factor(pdat$Distribution,levels=c("Breakpoint number","Copy-number","CN changepoint","Breakpoints per chr arm","Oscilating CN length","Segment size"))

## plot all signature weight
pdf(file = "CerCNsig_component_weights.pdf", width = 6, height = 15, onefile = F)
ggplot(pdat,
       aes(x=interaction(Feature,Distribution),
           y=value,fill=Distribution,group=Distribution))+
  geom_col(position="dodge")+
  scale_x_discrete(labels=c(1:5,1:7,1:5,1:7,1:2,1:10))+
  theme(legend.position="bottom",
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  ylab("")+
  coord_cartesian(ylim=c(0,1))+
  ggtitle("HGSC Cervical")+
  xlab("Component number")+
  facet_wrap(~sig, ncol = 1)
dev.off()
