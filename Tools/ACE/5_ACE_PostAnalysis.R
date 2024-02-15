#!/usr/bin/env Rscript

library(base)
library(utils)
library(stringr)
library(ACE)
library(Biobase)
library(QDNAseq)
library(ggplot2)

# A model which have been revised together with Rascal and ichorCNA
models <- read.csv("/home/minerva/WGS/ACE/Data/MaNiLa/models.tsv", header=T, sep="\t", stringsAsFactors = F)


temp <- list.files(getwd(), pattern="*rds")

	
	object <- readRDS(temp)
	id <- str_sub(temp, end = -8)
	name <- models[models$sample %in% id, ]$name
	c <- models[models$sample %in% id, ]$cellularity
	p <- models[models$sample %in% id, ]$ploidy


	template <- objectsampletotemplate(object)
	ACEcall(template, title = paste(name, id), cellularity = c, ploidy = p, standard = 1, cap = 18)
	
## Modelfile contains at least two columns: the first specifying the sample names and the second specifying the cellularity. The third column is the ploidy of the samples. When omitted, it is assumed to be 2. The fourth column is the standard of the samples.
	
	postanalysisloop(object, "/home/minerva/WGS/ACE/Data/MaNiLa/models.tsv", outputdir = "./postanalysis", imagetype = 'png') 




