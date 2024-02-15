#!/usr/bin/env bash

conda init
conda activate sWGS1.0


#Download fastQ files using sra-tools
prefetch SRR3412599
fastq-dump -X 10 -Z SRR3412599 -O /home/minerva/WGS/Data #extract 10 raw reads in standard output (pairend shows in one long read) from low-pass WGS OC primary tumor 
fastq-dump --split-e SRR3412599 -O /home/minerva/WGS/Data
head SRR3412599_1.fastq SRR3412599_2.fastq
