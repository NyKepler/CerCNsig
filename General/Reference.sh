#!/usr/bin/env bash

conda init
conda activate sWGS1.0


#Download reference and create Index
mkdir /home/minerva/WGS/Reference
#Download reference and create Index (hg19 UCSC was used in previous paper)
wget -N ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz -O /home/minerva/WGS/Reference/bwaGRCh37.fa.gz #recommended by BWA and best for QDNASeq


#instal pv to show progress of gunzip
sudo apt-get install pv #'brew install pv' in OS 
pv /home/minerva/WGS/Reference/bwaGRCh37.fa.gz | gunzip /home/minerva/WGS/Reference/bwaGRCh37.fa

bwa index -p /home/minerva/WGS/Reference/bwaGRCh37 -a bwtsw /home/minerva/WGS/Reference/bwaGRCh37.fa

picard CreateSequenceDictionary -R /home/minerva/WGS/Reference/bwaGRCh37.fa -O /home/minerva/WGS/Reference/bwaGRCh37.dict #For GATK tool kits

samtools faidx /home/minerva/WGS/Reference/bwaGRCh37.fa #For GATK tool kits

#Download Known-Sites files contains SNP, SNV-Indels and No-Match
wget -N ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/common_all_20180423.vcf.gz -O /home/minerva/WGS/KnownSites/ncbiGRCh37.vcf.gz # NCBI for GATK hg19/37 B151 dbSNP

wget -N https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz -O /home/minerva/WGS/KnownSites/gatkGRCh37.vcf.gz # GATK hg19/37 dbSNP

wget -N https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0/Mills_and_1000G_gold_standard.indels.b37.vcf.gz -O /home/minerva/WGS/KnownSites/gatkGRCh37.indel.vcf.gz # GATK hg19/37 Indels

#Unzip Known-Sites files
pv /home/minerva/WGS/KnownSites/ncbiGRCh37.vcf.gz | gunzip /home/minerva/WGS/KnownSites/ncbiGRCh37.vcf

pv /home/minerva/WGS/KnownSites/gatkGRCh37.vcf.gz | gunzip /home/minerva/WGS/KnownSites/gatkGRCh37.vcf

pv /home/minerva/WGS/KnownSites/gatkGRCh37.indel.vcf.gz | gunzip /home/minerva/WGS/KnownSites/gatkGRCh37.indel.vcf

#Index feature files for GATK tool kits
gatk IndexFeatureFile -I ncbiGRCh37.vcf 
gatk IndexFeatureFile -I gatkGRCh37.vcf
gatk IndexFeatureFile -I gatkGRCh37.indel.vcf



