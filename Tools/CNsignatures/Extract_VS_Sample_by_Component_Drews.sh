#!/usr/bin/env bash

# conda activate R
# execute by ./home/researcher/CerCNsig/Script/Extract_VS_Sample_by_Components_Drews.sh 

dir=/home/researcher/CerCNsig
version=HGSC_VS_good
method=rascal
input=$dir/input/$version.csv

#Create sample info list with header and fill in location of each file into config file
mkdir $dir/Randomforest_Drews
trainlist=$dir/Randomforest_Drews/Sample.List.$method.$version.txt
printf "Patient\tSample\tLibrary\tSegfile\tCell\n" > $trainlist  

#Create Batch job
while IFS=, read -r Patient Group Subgroup Interval ID Type Sample Library Fastq Binsize Ploidy Cellularity Segfile
do

#' segment files from ACE+Rascal+ichorCNA
#path=$dir/absegments
#segfile="$Library"_absegments.csv

path=$dir/abs_segfile
segfile="$Segfile"_seg.tsv
 
  cd $dir/Randomforest_Drews
  cp $path/$segfile .  #Copy selected type of samples 
  cut -f 3,4,5,6 $segfile > "$Library"_seg.tsv #select only "Chromosome", "Start", "End" and "Segment_Mean2" in rascal result
  sed -i '1s/.*/chromosome\tstart\tend\tsegVal/' "$Library"_seg.tsv #rename column names
  printf "$Patient\t$Sample\t$Library\t$segfile\t$Cellularity\n" >> $trainlist 
  

done < <(tail -n +2 $input)

#Run Extract_VS_Sample_by_Component Rscript
cd $dir/Randomforest_Drews
Rscript /home/researcher/CerCNsig/Script/Extract_VS_Sample_by_Component_Drews.R









