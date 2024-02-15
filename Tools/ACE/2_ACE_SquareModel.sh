#!/usr/bin/env bash

#conda activate R
#Execue by ./WGS/ACE/Script/2_ACE_SquareModel.sh

seq="CTG_2023_209"
project="MaNiLa_CS3"
input=MaNiLa_CS3_input.csv
group=MaNiLa_CS3_VS

#Create a csv file for ACE squaremodel result of all samples 
bestfit=/home/minerva/WGS/ACE/Data/$project/bestfits.$group.tsv
printf "Library\tPloidy\tCellularity\tError\tStandard\tPenalty\tPenploidy\tMethod\tSample\n" > $bestfit

negcall=/home/minerva/WGS/ACE/Data/$project/negcalls.$group.tsv
printf "Library\tPloidy\tCellularity\tError\tStandard\tPenalty\tPenploidy\tMethod\tSample\n" > $negcall

poscall=/home/minerva/WGS/ACE/Data/$project/poscalls.$group.tsv
printf "Library\tPloidy\tCellularity\tError\tStandard\tPenalty\tPenploidy\tMethod\tSample\n" > $poscall

while IFS=, read -r patient type sample library fastq bin
do

pwd=/home/minerva/WGS/Data/$project/$patient
mkdir -p $pwd

## remove old analysis data rm -rf $pwd/ACE*

#Create batch job for different binsizes
swd=$pwd/ACE/$sample
mkdir -p $swd
cd $swd

## Copy rds files from all samples to working directory
cp $pwd/QDNAseq/$sample/$library.$bin"kb".sg.rds .

## Run ACE Rscript
mkdir -p ./squaremodel
Rscript /home/minerva/WGS/ACE/Script/1_ACE_SquareModel.R
rm *rds


version=v2 # check 1_ACE_SquareModel.R 
## Extract the lowest penalty with minimum error from individual files from all samples to one file
mv ./squaremodel/$library.$bin"kb".bestfits.tsv ./squaremodel/$library.$bin"kb".bestfits.$version.tsv
sed -n '2{p;q}' ./squaremodel/$library.$bin"kb".bestfits.$version.tsv > tmp.tsv
cat tmp.tsv | sed s/$/"\t"$sample/g >> $bestfit
rm tmp.tsv
## Extract all negative call (cellularity = 1) from individual files from all samples to one file
mv ./squaremodel/$library.$bin"kb".negcalls.tsv ./squaremodel/$library.$bin"kb".negcalls.$version.tsv
sed '1d' ./squaremodel/$library.$bin"kb".negcalls.$version.tsv > tmp.tsv
cat tmp.tsv | sed s/$/"\t"$sample/g >> $negcall
rm tmp.tsv
## Extract all positive call from individual files from all samples to one file
mv ./squaremodel/$library.$bin"kb".poscalls.tsv ./squaremodel/$library.$bin"kb".poscalls.$version.tsv
sed '1d' ./squaremodel/$library.$bin"kb".poscalls.$version.tsv > tmp.tsv
cat tmp.tsv | sed s/$/"\t"$sample/g >> $poscall
rm tmp.tsv


done < <(tail -n +2 /home/minerva/WGS/QDNAseq/Data/MaNiLa/$input)


