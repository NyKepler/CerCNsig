#!/usr/bin/env Rscript

library(base)
library(utils)
library(stringr)
library(R.cache)
library(digest)
library(ACE)
library(Biobase)
library(QDNAseq)
library(ggplot2)
future::plan("multisession", workers=8) #Multiprocess with 4 cores
options("QDNAseq::verbose" = NA) #Send verbose output to stdout instead of stderr


#Process copynumbersegmented data rds* from QDNAseq_PE analysis

#runACE(outputdir = "./CNestimation", penalty = 0.5, filetype='rds', ploidies = c(2,3,4), cap = 24, savereadcounts = TRUE, imagetype='png', printsummaries = 2)


temp <- list.files(getwd(), pattern="*rds")

	id <- str_sub(temp, end = -8)
	object <- readRDS(temp)
	penalties <- seq(0, 0.5, 0.01) # Penalty is range from 0 - 1 but we have validated that should not more than 0.5
	penploidies <- 0.5
	methods <- c("MAE") # use MAE only as we have already compared "RMSE"
	sample <- c()
	ploidy <- c()
	cellularity <- c()
	error <- c()
	standard <- c()
	penalty <- c()
	penploidy <- c()
	method <- c()
	
## validate sensitivity and specificity of the model.
for (a in penalties) {
  for (b in penploidies) {
    for (c in methods) {
    	
# Version 1 (tumor) :sm <- squaremodel(object, QDNAseqobjectsample = TRUE, ptop = 4.3, pbottom = 1.8, prows = 250, penalty = a, penploidy = b, cellularities = seq(1,100), method = c)
# Version 2 (Blood, Plasma, Endome and VS): 
sm <- squaremodel(object, QDNAseqobjectsample = TRUE, ptop = 5, pbottom = 1, prows = 100, penalty = a, penploidy = b, cellularities = seq(5,100), method = c) 
	
	sample <- append(sample, id)
	ploidy <- append(ploidy, sm$minimadf[1,1])
	cellularity <- append(cellularity, sm$minimadf[1,2])
	error <- append(error, sm$minimadf[1,3])
	standard <- append(standard, sm$standard)
	penalty <- append(penalty, sm$penalty)
	penploidy <- append(penploidy, sm$penploidy)
	method <- append(method, sm$method)

    }
  }
} 
	bestfits <- data.frame(sample, ploidy, cellularity, error, standard, penalty, penploidy, method)
	bestfits <- bestfits[order(bestfits$error),]
	write.table(bestfits, file.path(paste0("./squaremodel/", id, ".bestfits.tsv")), quote = FALSE, sep = "\t", na = "", row.names = FALSE)
	negcalls <-  bestfits[bestfits$cellularity == 1,]
	write.table(negcalls, file.path(paste0("./squaremodel/", id, ".negcalls.tsv")), quote = FALSE, sep = "\t", na = "", row.names = FALSE)
	poscalls <-  bestfits[bestfits$cellularity != 1,]
	write.table(poscalls, file.path(paste0("./squaremodel/", id, ".poscalls.tsv")), quote = FALSE, sep = "\t", na = "", row.names = FALSE)
	
	
	
