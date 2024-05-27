#!/usr/bin/env bash

#conda activate R
#execute by ./WGS/CNsignatures/Script/Pat_Signatures_ACE_1.sh 

dir=/home/minerva/WGS/CNsignatures/Data/MaNiLa
method=final
version=revised.2
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

