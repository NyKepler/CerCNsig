#!/usr/bin/env bash

# conda activate sWGS1.0
# Execue script ./WGS/ichorCNA/Script/1_HMMcopy_Batch.sh

seq="CTG_2023_209"
project="MaNiLa_CS3"
input=MaNiLa_CS3_input.csv
group=MaNiLa_CS3_VS

while IFS=, read -r patient type sample library fastq bin
do

pwd=/home/minerva/WGS/Data/$project/$patient
mkdir -p $pwd

#Create batch job for different binsizes
assay=$pwd/ichorCNA/$sample
mkdir -p $assay
cd $assay

## Copy bam files to working directory
cp /home/minerva/WGS/fastq/$project/$patient/$library/"$library"_BQSR.bam "$library".bam
cp /home/minerva/WGS/fastq/$project/$patient/$library/"$library"_BQSR.bam.bai "$library".bam.bai

#You should not use any bin size smaller than 500kb if you do not have a matched normal sample. If you use smaller bins, then germline CNVs can confound the tumor fraction estimation. It is usually better to simply use 1Mb bins because it offers the cleanest normalized read coverage signals for estimating the tumor fraction.

#readCounter --window 10000 --quality 20 --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" ./"$library".bam > ./$library.10kb.wig

readCounter --window 500000 --quality 20 --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" ./"$library".bam > ./$library.500kb.wig

readCounter --window 1000000 --quality 20 --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" ./"$library".bam > ./$library.1000kb.wig

rm *bam*

done < <(tail -n +2 /home/minerva/WGS/QDNAseq/Data/MaNiLa/$input)

