#!/usr/bin/env bash

sudo ls

patientID="MaNiLa_BRCA_115"


# cd /home/minerva/WGS/Data/MaNiLa_CS1/"$patientID"_QDNAseq_v0
cd /mnt/DATA/WGS/MaNiLa_CS1/BAM/"$patientID"_BWA
# cd /home/minerva/WGS/Data/MaNiLa_CS1/"$patientID"_BWA
for i in {33,34}
do
sampleID="CS1_"$i"_S"$i
cd $sampleID

#Remove files before backup
# rm -f $sampleID.* "$sampleID"_marked*
# rm -f "$sampleID"_BQSR* "$sampleID"_dup*
# sudo rm -f "$sampleID".bam "$sampleID".sam "$sampleID".sorted.sam "$sampleID".sorted.mapped.bam

#Copy files for re-analysis
mkdir /home/minerva/WGS/Data/MaNiLa_CS1/"$patientID"_BWA/
mkdir /home/minerva/WGS/Data/MaNiLa_CS1/"$patientID"_BWA/$sampleID
sudo cp "$sampleID"_BQSR* "$sampleID"_dup* /home/minerva/WGS/Data/MaNiLa_CS1/"$patientID"_BWA/$sampleID


cd ..

echo $i
done

sudo chmod ugo+rwx /home/minerva/WGS/Data/MaNiLa_CS1/"$patientID"_BWA #Give full permission

# sudo cp -R /home/minerva/WGS/Data/MaNiLa_CS1/"$patientID"_QDNAseq_v0 /mnt/DATA/CTG_2021_047/MaNiLa_CS1/CNA

# sudo cp -R /home/minerva/WGS/Data/MaNiLa_CS1/"$patientID"_BWA /mnt/DATA/CTG_2021_047/MaNiLa_CS1/BAM

# sudo mv MaNiLa_HGSOC_140 MaNiLa_HGSOC_140_binsize_XY #change folder name
