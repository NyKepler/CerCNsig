#!/usr/bin/env bash

# conda activate R
# execute by ./WGS/CNsignatures/Script/Signatures_ACE.sh 

dir=/home/minerva/WGS/CNsignatures/Data/MaNiLa
method=final
version=revised.2
input=$dir/input.$method.$version.csv

#Create sample info list with header and fill in location of each file into config file
mkdir $dir/Tumor
tlist=$dir/Tumor/Tumor.List.$method.$version.txt
printf "Patient\tSample\tLibrary\tSegfile\tCell\n" > $tlist  

mkdir $dir/Cytology
clist=$dir/Cytology/Cytology.List.$method.$version.txt
printf "Patient\tSample\tLibrary\tSegfile\tCell\n" > $clist 

#Create Batch job
while IFS=, read -r Project Patient Type Sample Library Ploidy Cellularity
do
project=$Project
patient=$Patient
type=$Type
sample=$Sample
library=$Library
ploidy="$Ploidy"
cell=$Cellularity
path=/home/minerva/WGS/Data/$project/$patient/ACE/$sample/postanalysis/segmentfiles #change path for different assay
segfile="$library"_segments.tsv
#cna=`awk -F '\t' '{print $8}' $path/$segfile | tail -n +2 | sort | uniq -c | awk '$2 !=$ploidy {sum+=$1} END{print sum}'`
#cn=`awk -F '\t' '{print $8}' $path/$segfile | tail -n +2 | sort | uniq -c | awk '{sum+=$1} END{print sum}'`
#cin="$((100 * cna / cn))"


if [[ $patient == *"HGSOC"* ]] && [[ $type == *"Tumor"* ]] #&& [[ $cell -ge "10" ]] #|| $type == *"Tissue"* ]] #&& [[ $cell != "0" ]]
then
  cd $dir/Tumor
  cp $path/$segfile .  #Copy selected type of samples
  cut -f 1,2,3,6 $segfile > temp.tsv #select only "Chromosome", "Start", "End" and "Segment_Mean2" in ACE result
  sed -i '1s/.*/chromosome\tstart\tend\tsegVal/' temp.tsv #rename column names
  mv temp.tsv $segfile
  printf "$patient\t$sample\t$library\t$segfile\t$cell\n" >> $tlist
 
 
elif [[ $patient == *"HGSOC"* ]] && [[ $type == *"Tissue"* || $type == *"VS"* || $type == *"Endo"* ]] #&& [[ $cell != "0" ]] #Samples with either at least 20 CNA events or 1/3 segments have CNA are selected  
then 
  cd $dir/Cytology
  cp $path/$segfile .  #Copy selected type of samples 
  cut -f 1,2,3,6 $segfile > temp.tsv #select only "Chromosome", "Start", "End" and "Segment_Mean2" in ACE result
  sed -i '1s/.*/chromosome\tstart\tend\tsegVal/' temp.tsv #rename column names
  mv temp.tsv $segfile
  printf "$patient\t$sample\t$library\t$segfile\t$cell\n" >> $clist 
  
fi

done < <(tail -n +2 $input)

#Run CNsignature Rscript

cd $dir/Tumor
#Rscript /home/minerva/WGS/CNsignatures/Script/Extract_VS_signatures.R
#Rscript /home/minerva/WGS/CNsignatures/Script/Extract_tumor_signatures.R
Rscript /home/minerva/WGS/CNsignatures/Script/CNsignatures.R 
Rscript /home/minerva/WGS/CNsignatures/Script/Cytosignatures.R 

#cd $dir/Cytology
#Rscript /home/minerva/WGS/CNsignatures/Script/Extract_VS_signatures.R
#Rscript /home/minerva/WGS/CNsignatures/Script/Cytosignatures.R 
#Rscript /home/minerva/WGS/CNsignatures/Script/CNsignatures.R 
#Rscript /home/minerva/WGS/CNsignatures/Script/Extract_tumor_signatures.R





