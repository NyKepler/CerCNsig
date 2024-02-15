#!/usr/bin/env bash

## source /home/minerva/miniconda3/bin/activate
## conda init
## conda activate sWGS1.0 
## bash ./WGS/QDNAseq/Script/Alignment_1kg_h19.sh



pwd=/mnt/DATA/WGS/1000genomes/fastq
input=$pwd/input.csv

while IFS=, read -r Sample_Name Run_ID
do

case=$Sample_Name
run=$Run_ID

# Using sample id to extract fastq file name in sample_file_map.txt
  threads=8
  ref=/home/minerva/WGS/QDNAseq/Reference/bwaGRCh37 #NCBI bwa-indexed hg19 reference
  fastq=$pwd/fastq/$run*.filt.fastq.gz
  bam=$pwd/BAM/$run.sort.bam
  
if [[ ! -f $bam ]]
then
    
  bwa mem -t $threads -M -R '@RG\tID:"$case"\tSM:"$case"\tPL:ILLUMINA\tLB:Library' $ref $fastq | samtools sort -@ $threads > $bam
  
  samtools index $bam

#Extract comprehensive statistics of orginal bam file
  samtools stats -@ $threads $bam > /$pwd/BAM/$run.sort.stats.txt	 
  

fi

done < <(tail -n +2 $input)



