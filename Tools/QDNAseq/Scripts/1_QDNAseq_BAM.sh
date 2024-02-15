#!/usr/bin/env bash
# source /home/minerva/miniconda3/bin/activate
# conda init
# conda activate sWGS1.0 
# bash ./WGS/QDNAseq/Script/1_QDNAseq_BAM.sh

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
swd=$pwd/$library
mkdir -p $swd

ref=/home/minerva/WGS/QDNAseq/Reference/bwaGRCh37

bwa mem -t $threads -M -R '@RG\tID:"$library"\tSM:"$library"\tPL:ILLUMINA\tLB:Library' $ref $fastq | samtools sort -@ $threads > $swd/$library.sort.bam

rm -f $fastq

#Generate BAM index
samtools index $swd/$library.sort.bam

#samtools view -H $swd/"$library".sort.bam to check header

#Extract comprehensive statistics of orginal bam file
samtools stats -@ $threads $swd/$library.sort.bam > $swd/$library.sort.stats.txt

#Simple Statistics using samtools flagstat
samtools flagstat -@ $threads $swd/$library.sort.bam > $swd/$library.sort.fstat.txt

#Estimate library complexity from the sequence of read pairs (different to Picard)
picard EstimateLibraryComplexity -I $swd/"$library".sort.bam -O $swd/"$library"_complex_metrics.txt --MAX_RECORDS_IN_RAM 50000000 #no output file indicates complexity cannot be estimated

#Remove duplicate using Picard
gatk MarkDuplicatesSpark -I $swd/"$library".sort.bam -O $swd/"$library"_marked_dup.bam -M $swd/"$library"_dup_metrics.txt --remove-all-duplicates --spark-master local[$threads]

#Base Quality Score Recalibration using GATK
mkdir -p $swd/tmp
gatk BQSRPipelineSpark -R /home/minerva/WGS/QDNAseq/Reference/bwaGRCh37.fa -I $swd/"$library"_marked_dup.bam --known-sites /home/minerva/WGS/QDNAseq/KnownSites/gatkGRCh37.vcf --known-sites /home/minerva/WGS/QDNAseq/KnownSites/ncbiGRCh37.vcf --known-sites /home/minerva/WGS/QDNAseq/KnownSites/gatkGRCh37.indel.vcf -O $swd/"$library"_BQSR.bam --tmp-dir $swd/tmp --spark-master local[$threads]

#Index BAM file
samtools index $swd/"$library"_BQSR.bam


#remove tmp folder
rm -rf $swd/tmp

done < <(tail -n +2 /home/minerva/WGS/fastq/$project/$input)
