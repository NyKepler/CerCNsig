#!/usr/bin/env Rscript

library(base)
library(utils)
library(stringr)
library(R.cache)
library(digest)
library(Biobase)
library(ggplot2)
library(rascal)
library(dplyr)
future::plan("multisession", workers=4) #Multiprocess with 4 cores

dir="/home/minerva/WGS/Rascal/Data/MaNiLa/"

## Process copy number data extracted from QDNAseq_PE analysis
input <- read.csv(paste0(dir, "MaNiLa_CS3_input_rerun.csv"), header=T, sep=",", stringsAsFactors = F)

temp <- list.files(getwd(), pattern="*csv")

	name <- str_sub(temp, end = -5)
	sample <- input[input$Library %in% name, ]$Sample
	copy_number <- read.csv(temp, header=T, sep=",", stringsAsFactors = F)
	copy_number <-  mutate(copy_number, chromosome=factor(chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")))
	segments <- copy_number_segments(copy_number)
	
	write.table(segments, file.path(paste0("./segmentfiles/", name, "_segments.csv")), quote = FALSE, sep = ",", na = "", row.names = FALSE)
	
	
## Import solution result from Rascal 
solutions <- read.csv(paste0(dir, "Solutions.csv"), header = T, sep = ",", stringsAsFactors = F)

	p <- solutions[solutions$Library %in% name, ]$Ploidy
	c <- solutions[solutions$Library %in% name, ]$Cellularity
	d <- solutions[solutions$Library %in% name, ]$Distance


## Visualizing copy number profiles
my_theme<-ggplot2::theme_bw()+ggplot2::theme(axis.text=ggplot2::element_text(size=10),
                   axis.title=ggplot2::element_text(size=10),
                   strip.text.x = ggplot2::element_text(size = 14),
                   strip.text.y = ggplot2::element_text(size = 14),
                   legend.text = ggplot2::element_text(size = 14),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())

title <- paste(name, sample)               
	
	png(file = paste0(name, ".cn.png"), width = 800, height = 633)
	genome_copy_number_plot(copy_number, segments, sample = name, ylabel = "relative copy number") + my_theme + ggtitle(title)
	dev.off()
	
	#chromosome_copy_number_plot(copy_number, segments, sample = name, chromosome = 3)
	#chromosome_copy_number_plot(copy_number, segments, sample = name, chromosome = 3, start = 100000000, end = 120000000)
	
	# Copy number log2 ratio
	log2_ratio <- mutate(copy_number, copy_number = log2(copy_number))
	log2_ratio_segments <- mutate(segments, copy_number = log2(copy_number))
	write.table(log2_ratio_segments, file.path(paste0("./segmentfiles/", name, "_log2segments.csv")), quote = FALSE, sep = ",", na = "", row.names = FALSE)
	
	png(file = paste0(name, ".log2cn.png"), width = 800, height = 633)
	genome_copy_number_plot(log2_ratio, log2_ratio_segments, sample = name, min_copy_number = -2, max_copy_number = 3, xlabel = NULL, ylabel = expression(log[2]~ratio)) + my_theme + ggtitle(title)
	dev.off()
	
	#chromosome_copy_number_plot(log2_ratio, log2_ratio_segments, sample = name, chromosome = 3, xlabel = NULL, ylabel = expression(log[2]~ratio))
	#chromosome_copy_number_plot(copy_number, segments, sample = name, chromosome = 3, start = 100000000, end = 120000000)
	
	# Absolute copy number
	absolute_copy_number <- mutate(copy_number, across(c(copy_number, segmented), relative_to_absolute_copy_number, p, c))
	absolute_segments <- mutate(segments, copy_number = relative_to_absolute_copy_number(copy_number, p, c))
	write.table(absolute_segments, file.path(paste0("./segmentfiles/", name, "_absegments.csv")), quote = FALSE, sep = ",", na = "", row.names = FALSE)
	
	png(file = paste0(name, ".abcn.png"), width = 800, height = 633)
	genome_copy_number_plot(absolute_copy_number, absolute_segments, min_copy_number = 0, max_copy_number = 7, copy_number_breaks = 0:7, point_colour = "grey40", ylabel = "absolute copy number") + my_theme + ggtitle(title)
	dev.off()
	
	#chromosome_copy_number_plot(absolute_copy_number, absolute_segments, chromosome = 3, copy_number_breaks = 0:8, ylabel = "absolute copy number")
	

	png(file = paste0(name, ".abdensity.png"), width = 800, height = 633)
	copy_number_density_plot(absolute_copy_number$segmented, min_copy_number = 0, max_copy_number = 6, xlabel = "scaled/absolute copy number") + ggtitle(title)
	dev.off()
	
	copy_number_steps <- tibble(absolute_copy_number = 0:5)
	copy_number_steps <- mutate(copy_number_steps, copy_number = absolute_to_relative_copy_number(absolute_copy_number, p, c))
	
	png(file = paste0(name, ".density.png"), width = 800, height = 633)
	copy_number_density_plot(copy_number$segmented, copy_number_steps, min_copy_number = 0.3, max_copy_number = 1.7, xlabel = "relative copy number") + ggtitle(title)
	dev.off()
	
	
## The function `absolute_copy_number_distance_heatmap` takes relative copy numbers and computes the distance grid and generates the heatmap in one step. The `absolute_copy_number_distance` function can compute either a mean absolute deviation/distance (MAD) or root mean square deviation/distance (RMSD). 
	
	png(file = paste0(name, ".MAD.dishm.png"), width = 800, height = 633)
	absolute_copy_number_distance_heatmap(segments$copy_number, segments$weight, distance_function = "MAD") + ggtitle(title)
	dev.off()
	
## Use estimate TP53 MAF to adjust solutions which works better at smaller binsize like 30kb
	#TP53_cn <- filter(copy_number, chromosome == "17", end > 7565097, start < 7590856, segments != "NA")
	#TP53_solutions <- solutions %>% select(ploidy, cellularity) %>% mutate(tp53_absolute_copy_number = relative_to_absolute_copy_number(0.832, ploidy, cellularity)) %>% mutate(tp53_tumour_fraction = tumour_fraction(tp53_absolute_copy_number, cellularity))
	#write.table(TP53_solutions, file.path(paste0("./Solutions/", name, ".TP53.Solutions.csv")), quote = FALSE, sep = ",", na = "", row.names = FALSE)

