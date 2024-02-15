#!/usr/bin/env bash

#This is the 1st step of the BICseq workflow
#conda activate sWGS1.0
#exe by ./WGS/BICseq2/Script/1_BICseq2_Alignment.sh

seq="CTG_2023_209"
project="MaNiLa_CS3"
input=MaNiLa_CS3_input.csv

while IFS=, read -r patient type sample library fastq bin
do

pwd=/home/minerva/WGS/fastq/$project/$patient
mkdir -p $pwd

#Unzip fastq files from backup HDD disk, "sudo ls" to unlock permission before running
fastq1=/home/minerva/WGS/fastq/$seq/SeqOnly/fastq/HMLYHDMXY/$fastq*R1_001.fastq.gz
fastq2=/home/minerva/WGS/fastq/$seq/SeqOnly/fastq/HMLYHDMXY/$fastq*R2_001.fastq.gz
cat $fastq1 > $pwd/"$fastq"_R1_001.fastq.gz
cat $fastq2 > $pwd/"$fastq"_R2_001.fastq.gz
fastq=$pwd/"$fastq"*.fastq.gz
threads=32


#Alignment with BWA
swd=$pwd/$library/BICseq2
mkdir -p $swd

ref=/home/minerva/WGS/WGS_v2/BICseq2/Reference/hg19

#Alignment with BWA and export bam file (Output all found alignments for single-end or unpaired paired-end reads and flag with secondary alignment)
bwa mem -t $threads -a -R '@RG\tID:"$library"\tSM:"$library"\tPL:ILLUMINA\tLB:Library' $ref $fastq | samtools sort -@ $threads > $swd/$library.sort.bam

rm -f $fastq

#Generate BAM index
samtools index $swd/$library.sort.bam

#Extract comprehensive statistics of orginal bam file
samtools stats -@ $threads $swd/$library.sort.bam > $swd/$library.sort.stats.txt

#samtools view -H ./$sampleID/"$sampleID".sort.bam to check header

#Filter unique mapped read-pairs with same read length at least 100bp
samtools view -@ $threads -F 260 -f 3 -q 1 -h $swd/$library.sort.bam | awk '/^@/ || length($10) >= 100' | samtools view -@ $threads -Sb > $swd/$library.filt.bam

#Index filt bam file
samtools index $swd/$library.filt.bam

#Extract comprehensive statistics of filtered bam file
samtools stats -@ $threads $swd/$library.filt.bam > $swd/$library.filt.stats.txt

done < <(tail -n +2 /home/minerva/WGS/fastq/$project/$input)



#Run multiQC when all cases are aligned.
