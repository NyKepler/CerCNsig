#!/usr/bin/env bash

#conda activate R
#Execue by ./WGS/Rascal/Script/2_Rascal_complete.sh

seq="CTG_2023_209"
project="MaNiLa_CS3"
input=MaNiLa_CS3_input_rerun.csv # copy the input file from QDNAseq analysis and copy this name to the 1_Rascal_complete.R
group=MaNiLa_CS3_VS


#Create a txt file for Rascal solutions of all samples
solutions=/home/minerva/WGS/Rascal/Data/MaNiLa/Solutions.csv
printf "Library,Ploidy,Cellularity,Distance,Sample\n" > $solutions 

while IFS=, read -r patient type sample library fastq bin
do

pwd=/home/minerva/WGS/Data/$project/$patient
mkdir -p $pwd

#Create batch job for different binsizes
swd=$pwd/Rascal/$sample
mkdir -p $swd
cd $swd

## Copy rds files from all samples to working directory
cp $pwd/QDNAseq/$sample/$library.$bin"kb".sg.rds $library.rds


#Run Rascal Rscript to extract median-scaled copy number and solutions files. Distance model can be changed between MAD and RMSD.
mkdir -p solutions
version=default # check 1_Rascal_complete.R 
Rscript /home/minerva/WGS/Rascal/Script/extract_qdnaseq_copy_number_data.R -i $library.rds -m -o $library.csv
Rscript /home/minerva/WGS/Rascal/Script/fit_absolute_copy_number.R -i $library.rds -o ./solutions/$library.solutions.$version.txt --min-ploidy=2 --max-ploidy=5 --min-cellularity=0.05 
rm *rds

#extract solutions with lowest distance from individual files from all samples to one file
sed -n '2{p;q}' ./solutions/$library.solutions.$version.txt > temp.txt
cat temp.txt | sed s/$/","$sample/g >> $solutions
rm -f temp.txt

#Run Rascal Rscript to extract detail CN profiles and segment files
mkdir -p segmentfiles
Rscript /home/minerva/WGS/Rascal/Script/1_Rascal_complete.R 


#Move profile files to sample folder
mkdir -p profiles
mv *png profiles


done < <(tail -n +2 /home/minerva/WGS/Rascal/Data/MaNiLa/$input)

mv $solutions /home/minerva/WGS/Rascal/Data/MaNiLa/Solutions/Solutions.$group.$version.csv


