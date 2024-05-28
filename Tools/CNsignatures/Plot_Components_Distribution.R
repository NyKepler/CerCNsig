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
folder.name <- "Randomforest_CNsig_filt"
setwd(paste(path_home_r(), "CerCNsig", folder.name, sep = "/"))

##' import CN Components which extracted from HGSC_VS_good with CNsig over 3
CN_components <- readRDS("HGSC_VS_good_CNsig_over_3.Components.rds")

##' Distribution of 36 components that were extracted from the mixture modeling step
fillcol<- hue_pal()(6)
line_size <-0.7
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
dist_theme<-dist_theme+theme(axis.line = element_line(color = 'black',size=line_size))


#segsize
plotparam<-flexmix::parameters(CN_components[["segsize"]])
plotparam<-plotparam[,order(plotparam[1,])]
write.table(plotparam, "Components_segsize.csv", row.names = T, col.names = T)
ss<-ggplot(data = data.frame(x = c(1000,140000000)), aes(x)) + ylab("")+my_theme
for(i in 1:ncol(plotparam))
{
  ss<-ss+stat_function(fun = dnorm, n = 1000, args = list(mean = plotparam[1,i], sd = plotparam[2,i]),color=cbPalette[i],size=line_size)
}
ss<-ss+scale_x_continuous(breaks=c(30000000, 100000000))+xlab("Segment size")

ss

#breakpoint
plotparam<-flexmix::parameters(CN_components[["bp10MB"]])
bp<-ggplot(data = data.frame(x = c(0:10)), aes(x)) + ylab("")+theme()+scale_x_continuous(breaks=c(0,4,8))+my_theme+
  stat_function(geom="line",n=11,fun = dpois,args=list(lambda = plotparam[1]),color=cbPalette[1],size=line_size) +
  stat_function(geom="line",n=11,fun = dpois, args=list(lambda = plotparam[2]),color=cbPalette[2],size=line_size) +
  stat_function(geom="line",n=11,fun = dpois, args=list(lambda = plotparam[3]),color=cbPalette[3],size=line_size)+xlab("Breakpoint count per 10MB")
bp

plotparam<-plotparam[order(plotparam)] %>% as.data.frame()
write.table(plotparam, "Components_bp10MB.csv", row.names = T, col.names = T)

#oscilating
plotparam<-flexmix::parameters(CN_components[["osCN"]])
os<-ggplot(data = data.frame(x = c(0:10)), aes(x)) + ylab("")+
  scale_x_continuous(breaks=c(0,4,8))+my_theme+
  stat_function(geom="line",n=11,fun = dpois,args=list(lambda = plotparam[1]),color=cbPalette[1],size=line_size) +
  stat_function(geom="line",n=11,fun = dpois, args=list(lambda = plotparam[2]),color=cbPalette[2],size=line_size) +
  stat_function(geom="line",n=11,fun = dpois, args=list(lambda = plotparam[3]),color=cbPalette[3],size=line_size)+xlab("Length of oscilating CN segments")
os
plotparam<-plotparam[order(plotparam)] %>% as.data.frame()
write.table(plotparam, "Components_osCN.csv", row.names = T, col.names = T)


#changepoint
plotparam<-flexmix::parameters(CN_components[["changepoint"]])
plotparam<-plotparam[,order(plotparam[1,])]
write.table(plotparam, "Components_changepoint.csv", row.names = T, col.names = T)
cp<-ggplot(data = data.frame(x = c(0,3)), aes(x)) + ylab("")+
  scale_x_continuous(breaks=c(1,2,3))+my_theme
for(i in 1:ncol(plotparam))
{
  cp<-cp+stat_function(fun = dnorm, n = 1000, args = list(mean = plotparam[1,i], 
                                                          sd = plotparam[2,i]),color=cbPalette[i],size=line_size)
}
cp<-cp+xlab("Copy-number change point")


#copy-number
plotparam<-flexmix::parameters(CN_components[["copynumber"]])
plotparam<-plotparam[,order(plotparam[1,])]
write.table(plotparam, "Components_copynumber.csv", row.names = T, col.names = T)
cn<-ggplot(data = data.frame(x = c(1.5,3)), aes(x)) + ylab("")+
  scale_x_continuous(breaks=c(2,3))+my_theme
for(i in 1:ncol(plotparam))
{
  cn<-cn+stat_function(fun = dnorm, n = 1000, args = list(mean = plotparam[1,i],sd = plotparam[2,i]),color=cbPalette[i],size=line_size)
}
cn<-cn+xlab("Copy-number")


#bp per chr arm
plotparam<-flexmix::parameters(CN_components[["bpchrarm"]]) 
ct<-ggplot(data = data.frame(x = c(0:35)), aes(x)) + ylab("")+
  scale_x_continuous(breaks=c(0,10,30))+my_theme+
  stat_function(geom="line",n=36,fun = dpois,args=list(lambda = plotparam[1]),color=cbPalette[1],size=line_size) +
  stat_function(geom="line",n=36,fun = dpois, args=list(lambda = plotparam[2]),color=cbPalette[2],size=line_size) +
  stat_function(geom="line",n=36,fun = dpois, args=list(lambda = plotparam[3]),color=cbPalette[3],size=line_size) +
  stat_function(geom="line",n=36,fun = dpois, args=list(lambda = plotparam[4]),color=cbPalette[4],size=line_size) +
  stat_function(geom="line",n=36,fun = dpois, args=list(lambda = plotparam[5]),color=cbPalette[5],size=line_size)+xlab("Breakpoint count per chr arm") 
ct

plotparam<-plotparam[order(plotparam)] %>% as.data.frame()
write.table(plotparam, "Components_bpchrarm.csv", row.names = T, col.names = T)


pdf(file = paste0("HGSC_VS_filt_36_components_distribution.pdf"), width = 18, height = 3)
plot_grid(bp,cn,cp,ct,os,ss,ncol = 6)
dev.off()