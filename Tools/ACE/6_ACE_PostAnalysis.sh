#!/usr/bin/env bash

#conda activate R
#Execue by ./WGS/ACE/Script/ACE_PostAnalysis.sh

project="MaNiLa_CS3"
group=MaNiLa_CS3_VS
version=benchmark.revised
input=/home/minerva/WGS/ACE/Data/MaNiLa/input.$group.$version.csv 
#remember to fix the models.tsv file for the Rscript

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
rm -rf ./postanalysis*
mkdir -p ./postanalysis
Rscript /home/minerva/WGS/ACE/Script/5_ACE_PostAnalysis.R

rm *rds

#rename ACEcall plot and move into postanalysis folder
mv ./Rplots.pdf ./postanalysis/$library.$version.pdf 
convert -resize 1000x807 -quality 0 ./postanalysis/$library.$version.pdf ./postanalysis/$library.$version.png

done < <(tail -n +2 $input)




