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
future::plan("multisession", workers=4) #Multiprocess with 4 cores
options("QDNAseq::verbose" = NA) #Send verbose output to stdout instead of stderr

input <- read.csv("/home/minerva/WGS/ACE/Data/MaNiLa/input.MaNiLa_CS3_VS.lowest.penalty.csv", header=T, sep=",", stringsAsFactors = F)
temp <- list.files(getwd(), pattern="*rds")

	id <- str_sub(temp, end = -8)
	object <- readRDS(temp)
	penalty <- input[input$library %in% id, ]$penalty
	method <- input[input$library %in% id, ]$method
	model <- input[input$library %in% id, ]$model
	penploidy <- input[input$library %in% id, ]$penploidy 
		
if (model == "v1") {
	sm <- squaremodel(object, QDNAseqobjectsample = TRUE, method = method, ptop = 4.3, pbottom = 1.8, prows = 250, penalty = penalty, penploidy = penploidy, cellularities = seq(1,100))
} else {
	sm <- squaremodel(object, QDNAseqobjectsample = TRUE, method = method, ptop = 5, pbottom = 1, prows = 100, penalty = penalty, penploidy = penploidy, cellularities = seq(5,100)) 
}
	
	squaremodelsummary(object, QDNAseqobjectsample = TRUE, squaremodel = sm, samplename = id, printplots = TRUE, outputdir = "./squaremodel", imagetype = 'png')
	
	library <- id
	method <- sm$method
	penalty <- sm$penalty
	penploidy <- sm$penploidy
	standard <- sm$standard
	result <- head(sm$minimadf, 7)
	bestfits <- data.frame(library, method, penalty, penploidy, standard, result)
	write.table(bestfits, file.path(paste0("./squaremodel/", id, ".7bestfits.tsv")), quote = FALSE, sep = "\t", na = "", row.names = FALSE)


	
