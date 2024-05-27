##' Distribution of 39 components that were extracted from the mixture modeling step
dist_theme<-dist_theme+theme(axis.line = element_line(color = 'black',size=line_size))
load("~/CerCNsig/Randomforest/HGSC_Tumor.RData")
CN_components <- VS_components
line_size <-2

#segsize
plotparam<-flexmix::parameters(CN_components[["segsize"]])
plotparam<-plotparam[,order(plotparam[1,])]
ss<-ggplot(data = data.frame(x = c(1000,119118142)), aes(x)) + ylab("")+my_theme
for(i in 1:ncol(plotparam))
{
  ss<-ss+stat_function(fun = dnorm, n = 1000, args = list(mean = plotparam[1,i], sd = plotparam[2,i]),color=cbPalette[i],size=line_size)
}
ss<-ss+scale_x_continuous(breaks=c(100000000))+xlab("Segment size")

ss

#breakpoint
plotparam<-flexmix::parameters(CN_components[["bp10MB"]])
bp<-ggplot(data = data.frame(x = c(0:10)), aes(x)) + ylab("")+theme()+scale_x_continuous(breaks=c(0,4,8))+my_theme+
  stat_function(geom="line",n=11,fun = dpois,args=list(lambda = plotparam[1]),color=cbPalette[1],size=line_size) +
  stat_function(geom="line",n=11,fun = dpois, args=list(lambda = plotparam[2]),color=cbPalette[2],size=line_size) +
  stat_function(geom="line",n=11,fun = dpois, args=list(lambda = plotparam[3]),color=cbPalette[3],size=line_size)+xlab("Breakpoint count per 10MB")
bp

#oscilating
plotparam<-flexmix::parameters(CN_components[["osCN"]])
os<-ggplot(data = data.frame(x = c(0:10)), aes(x)) + ylab("")+
  scale_x_continuous(breaks=c(0,4,8))+my_theme+
  stat_function(geom="line",n=11,fun = dpois,args=list(lambda = plotparam[1]),color=cbPalette[1],size=line_size) +
  stat_function(geom="line",n=11,fun = dpois, args=list(lambda = plotparam[2]),color=cbPalette[2],size=line_size) +
  stat_function(geom="line",n=11,fun = dpois, args=list(lambda = plotparam[3]),color=cbPalette[3],size=line_size)+xlab("Length of oscilating CN segments")
os

#changepoint
plotparam<-flexmix::parameters(CN_components[["changepoint"]])
plotparam<-plotparam[,order(plotparam[1,])]
cp<-ggplot(data = data.frame(x = c(0,25)), aes(x)) + ylab("")+
  scale_x_continuous(breaks=c(4,12,22))+my_theme
for(i in 1:ncol(plotparam))
{
  cp<-cp+stat_function(fun = dnorm, n = 1000, args = list(mean = plotparam[1,i], 
                                                          sd = plotparam[2,i]),color=cbPalette[i],size=line_size)
}

cp<-cp+xlab("Copy-number change point")

#copy-number
plotparam<-flexmix::parameters(CN_components[["copynumber"]])
plotparam<-plotparam[,order(plotparam[1,])]
cn<-ggplot(data = data.frame(x = c(0,34)), aes(x)) + ylab("")+
  scale_x_continuous(breaks=c(4,16,32))+my_theme
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

pdf(file = paste0("HGSC_VS_39_components_distribution.pdf"), width = 18, height = 3)
plot_grid(bp,cn,cp,ct,os,ss,ncol = 6)
dev.off()
