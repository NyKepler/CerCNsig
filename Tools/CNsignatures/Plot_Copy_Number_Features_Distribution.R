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
library(data.table)


## 1. Import data
folder.name <- "Randomforest_CNsig_filt"
setwd(paste(path_home_r(), "CerCNsig", folder.name, sep = "/"))

##' supplementary figure to show distribution the copy number features
chrlen<-read.table("/home/researcher/CerCNsig/Script/data/hg19.chrom.sizes.txt",sep="\t",stringsAsFactors = F)[1:24,]
gaps<-read.table("/home/researcher/CerCNsig/Script/data/gap_hg19.txt",sep="\t",header=F,stringsAsFactors = F)

##' Prepare individual data files
load("All_Cervical.rdata")
All_Cervical <- read_csv("/home/researcher/CerCNsig/input/All_Cervical.csv")
HGSC.VS <- All_Cervical[All_Cervical$Group %like% "HGSC",] # 128 samples
Benign.VS <- All_Cervical[All_Cervical$Group %like% "Benign",] # 55 samples
RRSO.VS <- All_Cervical[All_Cervical$Group %like% "RRSO",] # 30 samples

segment_tables.HGSC.VS <- segment_tables[c(HGSC.VS$Library)]
segment_tables.Benign.VS <- segment_tables[c(Benign.VS$Library)]
segment_tables.RRSO.VS <- segment_tables[c(RRSO.VS$Library)]

load("HGSC_Tumor.rdata")
segment_tables.HGSC.Tumor <- segment_tables # 18 samples


#' set color and plot theme
fillcol<- hue_pal()(6)
line_size<-0.25
cbPalette <- c(RColorBrewer::brewer.pal(8,"Set2"),RColorBrewer::brewer.pal(9,"Set1"),"black")
ggplot2::theme_set(ggplot2::theme_gray(base_size = 10))
my_theme<-ggplot2::theme_bw()+ggplot2::theme(axis.text=ggplot2::element_text(size=10),
                                             axis.title=ggplot2::element_text(size=10),
                                             strip.text.x = ggplot2::element_text(size = 14),
                                             strip.text.y = ggplot2::element_text(size = 14),
                                             legend.text = ggplot2::element_text(size = 14),
                                             panel.grid.minor = ggplot2::element_blank(),
                                             panel.grid.major = ggplot2::element_blank())
dist_theme<-my_theme+
  theme(panel.border = element_blank(),
        plot.margin = unit(c(0,0.2,-0.3,-0.3), "cm"),
        axis.text=element_text(size=5),
        axis.ticks = element_line(colour = "black", size = 0.25),
        axis.ticks.length=unit(0.05, "cm"),
        plot.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(margin=ggplot2::margin(0,1,0,0,"pt")))

#' Extract segment features and distribution
segment_tables <- segment_tables.HGSC.Tumor # change for different sample type and group
sample.type <- "HGSC.Tumor" # change accordingly
all_dist_dat<-c()
seg_num<-getBPnum(segment_tables,chrlen)
all_dist_dat<-rbind(all_dist_dat,cbind(seg_num,distribution="BPnum"))

CN<-getCN(segment_tables)
all_dist_dat<-rbind(all_dist_dat,cbind(CN,distribution="copynum"))

cpCN<-getChangepointCN(segment_tables)
all_dist_dat<-rbind(all_dist_dat,cbind(cpCN,distribution="copynumcp"))

centromeres<-gaps[gaps[,8]=="centromere",]
centd<-getCentromereDistCounts(segment_tables,centromeres,chrlen)
colnames(centd)<-c("ID","value")
all_dist_dat<-rbind(all_dist_dat,cbind(centd,distribution="chrarm"))

segsize<-getSegsize(segment_tables)
all_dist_dat<-rbind(all_dist_dat,cbind(segsize,distribution="segsize"))

osc<-getOscilation(segment_tables)
all_dist_dat<-rbind(all_dist_dat,cbind(osc,distribution="osc"))

all_dist_dat$distribution<-plyr::revalue(all_dist_dat$distribution,
                                         c(BPnum="Breakpoint\ncount\nper 10MB",
                                           copynum="Copy\n number",
                                           copynumcp="CN\n changepoint",
                                           chrarm="Breakpoint\ncount\nper chr arm",
                                           segsize="Segment\n size",
                                           osc="Oscilating\n CN length"))

all_dist_dat$distribution=factor(all_dist_dat$distribution,levels = c("Breakpoint\ncount\nper 10MB","Copy\n number","CN\n changepoint","Breakpoint\ncount\nper chr arm","Oscilating\n CN length","Segment\n size"))

all_dist_dat<-data.frame(all_dist_dat,stringsAsFactors = F)

#plots for overview figures
bp<-ggplot(all_dist_dat[all_dist_dat$distribution=="Breakpoint\ncount\nper 10MB",],aes(x=as.numeric(value)))+
  geom_histogram(stat="count",size=line_size, fill=fillcol[1])+coord_cartesian(xlim=c(0,10))+xlab("")+ylab("")+dist_theme+
  scale_x_continuous(breaks=c(0,4,8))+ggtitle("Breakpoint count per 10MB")

cn<-ggplot(all_dist_dat[all_dist_dat$distribution=="Copy\n number",],aes(x=as.numeric(value)))+
  geom_density(size=line_size, fill=fillcol[2])+coord_cartesian(xlim=c(0,34))+xlab("")+ylab("")+dist_theme+
  scale_x_continuous(breaks=c(4,16,32))+ggtitle("Copy-number")

cp<-ggplot(all_dist_dat[all_dist_dat$distribution=="CN\n changepoint",],aes(x=as.numeric(value), fill=fillcol[5]))+
  geom_density(size=line_size, fill=fillcol[3])+coord_cartesian(xlim=c(0,25))+xlab("")+ylab("")+dist_theme+
  scale_x_continuous(breaks=c(4,12,22))+ggtitle("CN changepoint")

ct<-ggplot(all_dist_dat[all_dist_dat$distribution=="Breakpoint\ncount\nper chr arm",],aes(x=as.numeric(value)))+
  geom_histogram(stat="count",size=line_size, fill=fillcol[4])+coord_cartesian(xlim=c(0,35))+xlab("")+ylab("")+dist_theme+
  scale_x_continuous(breaks=c(0,10,30))+ggtitle("Breakpoint count per chr arm")

os<-ggplot(all_dist_dat[all_dist_dat$distribution=="Oscilating\n CN length",],aes(x=as.numeric(value)))+
  geom_histogram(stat="count",size=line_size, fill=fillcol[5])+coord_cartesian(xlim=c(0,10))+xlab("")+ylab("")+dist_theme+
  scale_x_continuous(breaks=c(0,4,8))+ggtitle("Oscilating CN chain length")

ss<-ggplot(all_dist_dat[all_dist_dat$distribution=="Segment\n size",],aes(x=as.numeric(value)))+
  geom_density(size=line_size, fill=fillcol[6])+dist_theme+coord_cartesian(xlim=c(0,119118142))+xlab("")+ylab("")+
  scale_x_continuous(breaks = c(100000000))+ggtitle("Segment size")

pdf(file = paste0(sample.type, ".CN.features.distribution.pdf"), width = 18, height = 3)
plot_grid(bp,cn,cp,ct,os,ss,ncol = 6)
dev.off()
