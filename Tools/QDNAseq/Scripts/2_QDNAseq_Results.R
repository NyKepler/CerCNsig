#!/usr/bin/env Rscript

library(base)
library(utils)
library(stringr)
library(Biobase)
library(QDNAseq)
source("/home/minerva/WGS/CNsignatures/Script/helper_functions.R")

#Choose the data file based on fitpicker from runACE.

temp <- list.files(getwd(), pattern="*rds")
name <- str_sub(temp, end = -5)
	
	object <- readRDS(temp)
	Library <- object@phenoData@data$name
	TotalRead <- object@phenoData@data$total.reads
	UsedRead <- object@phenoData@data$used.reads
	ExpVar <- object@phenoData@data$expected.variance
	Loess <- object@phenoData@data$loess.span
	filteredCN <- object[fData(object)$use,]
	segmented <- filteredCN[assayDataElement(filteredCN,"segmented")>0,]
	segtable <- getSegTable(segmented)
	Segments <- nrow(segtable)
	write.table(segtable, file.path(name, "/SegTable.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)
	
	Result <- data.frame(Library, TotalRead, UsedRead, Segments, ExpVar, Loess)
	write.table(Result, file.path(name, "/Results.csv"), quote = FALSE, sep = ",", na = "", row.names = FALSE)
	
	
	# Output a table of all segments with a probability of loss or gain is greater than 0.99
	filteredCN<-object[fData(object)$use,]
	regions_of_loss<-filteredCN[assayDataElement(filteredCN,"probloss") > 0.99,]
	relative_loss <- getSegTable(regions_of_loss)
	write.table(relative_loss, file = paste0(name, "/Relative_loss_p99.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
	regions_of_gain <- filteredCN[assayDataElement(filteredCN, "probgain") > 0.99, ]
	relative_gain <- getSegTable(regions_of_gain)
	write.table(relative_gain, file = paste0(name, "/Relative_gain_p99.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


