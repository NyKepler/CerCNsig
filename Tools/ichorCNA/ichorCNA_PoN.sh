#!/usr/bin/env bash

# conda activate R
# Execue script ./WGS/ichorCNA/Script/ichorCNA_PoN.sh


input=input_PoN.csv
pwd=/home/minerva/WGS/ichorCNA/Data/MaNiLa
mkdir -p $pwd
cd $pwd

while IFS=, read -r Project Patient Type Sample Library Fastq Bin
do

#Create our own Panel of Normal of individual sample types from the Benign group!   
project=$Project
patient=$Patient
type=$Type
sample=$Sample
library=$Library
bin=1000

if [[ "$type" == "ArchivalVS" ]]
then
wiglist=MaNiLa_Benign_"$type"_wig_files_"$bin"kb.txt
outfile=MaNiLa_Benign_"$type"_PoN_"$bin"kb
wig=/home/minerva/WGS/Data/$project/$patient/ichorCNA/$sample/$library."$bin"kb.wig
printf "$wig\n" >> $wiglist

fi

done < <(tail -n +2 $pwd/$input)

pkg=/home/minerva/miniconda3/envs/R/share/r-ichorcna-0.3.2-1

Rscript $pkg/scripts/createPanelOfNormals.R --filelist $pwd/$wiglist --gcWig $pkg/extdata/gc_hg19_"$bin"kb.wig --mapWig $pkg/extdata/map_hg19_"$bin"kb.wig --centromere $pkg/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt --outfile $outfile
