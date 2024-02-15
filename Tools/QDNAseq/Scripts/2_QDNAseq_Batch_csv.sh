#!/usr/bin/env bash

## conda activate R
## Execue script ./WGS/QDNAseq/Script/QDNAseq_Batch_csv.sh

seq="CTG_2023_209"
project="MaNiLa_CS3"
group=QDNAseq
input=MaNiLa_CS3_input_ichorCNA.csv # ! change the input file in the Rscript as well ! 

## Create a csv file for QDNAseq result of all samples

results=/home/minerva/WGS/QDNAseq/Data/MaNiLa/Results.$group.csv
printf "Library,TotalRead,UsedRead,Segments,ExpVar,Loess\n" > $results 

while IFS=, read -r patient type sample library fastq bin
do

pwd=/home/minerva/WGS/Data/$project/$patient
mkdir -p $pwd
#rm -rf $pwd/*PE* # remove old data if there is any

## Create batch job for different binsizes

assay=$pwd/QDNAseq/$sample
mkdir -p $assay
cd $assay

## Copy bam files to working directory
cp /home/minerva/WGS/fastq/$project/$patient/$library/"$library"_BQSR.bam "$library".bam
cp /home/minerva/WGS/fastq/$project/$patient/$library/"$library"_BQSR.bam.bai "$library".bam.bai


## Run QDNAseq copy number calling base on sample type
if [[ "$type" == *"Tumor"* ]]
then
 Rscript /home/minerva/WGS/QDNAseq/Script/QDNAseq_Tumor_PE_csv.R #remember to change the csv file path in the Rscript
else
 Rscript /home/minerva/WGS/QDNAseq/Script/QDNAseq_PE_csv.R #remember to change the csv file path in the Rscript
fi

rm *bam*

## Extract result from copy number data and save them to the 'Results' file
cd ..
cp $assay/$library.$bin"kb".cn.rds $sample.rds
Rscript /home/minerva/WGS/QDNAseq/Script/QDNAseq_Results.R
sed '1d' $assay/Results.csv >> $results
rm *rds

mv $assay/Results.csv $assay/Results.$bin"kb".csv
mv $assay/SegTable.tsv $assay/SegTable.$bin"kb".tsv
mv $assay/Relative_loss_p99.tsv $assay/Relative_loss_p99.$bin"kb".tsv
mv $assay/Relative_gain_p99.tsv $assay/Relative_gain_p99.$bin"kb".tsv


done < <(tail -n +2 /home/minerva/WGS/QDNAseq/Data/MaNiLa/$input)

mv $results /home/minerva/WGS/QDNAseq/Data/MaNiLa/Results/Results.$group.$bin"kb".csv






