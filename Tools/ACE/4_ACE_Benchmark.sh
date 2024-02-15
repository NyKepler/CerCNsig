#!/usr/bin/env bash

#conda activate R
#Execue by ./WGS/ACE/Script/4_ACE_Benchmark.sh

project="MaNiLa_CS3"
group=MaNiLa_CS3_VS
version=lowest.penalty
input=/home/minerva/WGS/ACE/Data/MaNiLa/input.$group.$version.csv # ! change the input file in the Rscript as well !

#Create a csv file for ACE squaremodel result of all samples
bestfit=/home/minerva/WGS/ACE/Data/MaNiLa/Bestfits.running.tsv
printf "Library\tMethod\tPenalty\tPenploidy\tStandard\tPloidy\tCellularity\tError\tMinimum\tSample\n" > $bestfit


while IFS=, read -r patient library penalty penploidy method sample model
do


pwd=/home/minerva/WGS/Data/$project/$patient
mkdir -p $pwd


#Create batch job for different binsizes
swd=$pwd/ACE/$sample
mkdir -p $swd
cd $swd

## Copy rds files from all samples to working directory
cp $pwd/QDNAseq/$sample/$library.sg.rds .

## Run ACE Rscript
mkdir -p ./squaremodel
Rscript /home/minerva/WGS/ACE/Script/3_ACE_Benchmark.R
rm *rds

## Extract the 7 cellularity models using the lowest penalty from individual samples to one file
mv ./squaremodel/$library.7bestfits.tsv ./squaremodel/$library.7bestfits.$version.tsv
sed -n '2{p;q}' ./squaremodel/$library.7bestfits.$version.tsv > tmp.tsv
cat tmp.tsv | sed s/$/"\t"$sample/g >> $bestfit
rm tmp.tsv

mv ./squaremodel/$library.png ./squaremodel/$library.$version.png
 
done < <(tail -n +2 $input)

## rename folder and files with binsize and version the squaremodel

mv $bestfit /home/minerva/WGS/ACE/Data/MaNiLa/Bestfits.$group.$version.tsv



