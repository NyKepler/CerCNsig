#!/usr/bin/env bash

# This is the 2nd step of the BICseq workflow
# conda activate sWGS1.0
# ./WGS/BICseq2/Script/2_BICseq_Preprocess.sh

#Define modified samtools helper script and run by perl
samtool_gu=/home/minerva/WGS/WGS_v2/BICseq2/modifiedSamtools/misc/samtools.pl

project="MaNiLa_CS3"
input=MaNiLa_CS3_input.csv

while IFS=, read -r patient type sample library fastq bin
do

pwd=/home/minerva/WGS/Data/$project/$patient/BICseq2/$sample/Seq_Out
mkdir -p $pwd


#Extract BWA uniquely-aligned reads with modified samtools and generate seq file

for x in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M}
do 
threads=32
readPosFile=$pwd/$library.chr$x.seq
bam=/home/minerva/WGS/fastq/$project/$patient/$library/BICseq2

#Generate seq file
samtools view -@ $threads $bam/$library.filt.bam chr$x | perl $samtool_gu unique - | cut -f 4 > $readPosFile

echo $x
done

done < <(tail -n +2 /home/minerva/WGS/fastq/$project/$input)
