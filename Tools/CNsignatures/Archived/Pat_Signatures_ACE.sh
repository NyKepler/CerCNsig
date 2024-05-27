#!/usr/bin/env bash

#conda activate R
#execute by ./WGS/CNsignatures/Script/Pat_Signatures_ACE.sh 

dir=/home/minerva/WGS/CNsignatures/Data/MaNiLa
method=final
version=revised
input1=$dir/input.$method.$version.csv
input2=$dir/input_pat_signature.csv

while IFS=, read -r Order Patient 
do
patient=$Patient
assay=$dir/$patient
mkdir -p $assay
cd $assay

#Create sample info list with header and fill in location of each file into config file
list=$patient.sample.list.txt
printf "Project\tPatient\tType\tSample\tLibrary\tSegfile\tCell\n" > $list

done < <(tail -n +2 $input2)

while IFS=, read -r Project Patient Type Sample Library Ploidy Cellularity
do
project=$Project
patient=$Patient
type=$Type
sample=$Sample
library=$Library
cell=$Cellularity
assay=$dir/$patient
cd $assay 

#Copy ACE segment files and prepare for CNsignature Rscript
path=/home/minerva/WGS/Data/$project/$patient/ACE/$sample/postanalysis/segmentfiles
segfile="$library"_segments.tsv
cp $segfile .  #Copy all samples from same patient
cut -f 1,2,3,6 "$library"_segments.tsv > temp.tsv #select only "Chromosome", "Start", "End" and "Segment_Mean2" in ACE result
sed -i '1s/.*/chromosome\tstart\tend\tsegVal/' temp.tsv #rename column names
mv temp.tsv $segfile
printf "$project\t$patient\t$type\t$sample\t$library\t$segfile\t$cell\n" >> $list
 
done < <(tail -n +2 $input1)


while IFS=, read -r Order Patient 
do
patient=$Patient
assay=$dir/$patient
cd $assay

#Run Cytosignature Rscript
Rscript /home/minerva/WGS/CNsignatures/Script/Cytosignatures.R 
Rscript /home/minerva/WGS/CNsignatures/Script/CNsignatures.R 

done < <(tail -n +2 $input2)
