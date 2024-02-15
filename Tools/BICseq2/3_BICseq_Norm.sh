#!/usr/bin/env bash

# This is the 3rd step of the BICseq workflow
# conda activate R
# execue by ./WGS/BICseq2/Script/3_BICseq_Norm.sh

project="MaNiLa_CS3"
input=MaNiLa_CS3_input.csv

while IFS=, read -r patient type sample library fastq bin
do

pwd=/home/minerva/WGS/Data/$project/$patient/BICseq2/$sample
mkdir -p $pwd

bin=10 #Default 100bp but according to the latest publication We normalized data and set the initial bin size as 10 bp for segmentation. The penalty parameter λ was chosen as 1.2
binsize=$bin"bp"
norm=$pwd/Norm_Out_$binsize
mkdir -p $norm


#Create config file with header and fill in location of each file into config file
config=$norm/$library.norm-config.$binsize.txt
printf "chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm\n" > $config  

for x in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M}
do 
chromName=chr$x
faFile=/home/minerva/WGS/WGS_v2/BICseq2/Reference/hg19.chr$x.fa
MapFile=/home/minerva/WGS/WGS_v2/BICseq2/MapFile/hg19.CRC.100mer.chr$x.txt
readPosFile=$pwd/Seq_Out/$library.chr$x.seq 
binFileNorm=$norm/$library.chr$x.$binsize.norm.bin
printf "$chromName\t$faFile\t$MapFile\t$readPosFile\t$binFileNorm\n" >> $config 

echo $x
done

#Normalizing potential biases in the sequencing data using BICseq2-norm
BICSEQ_NORM=/home/minerva/WGS/WGS_v2/BICseq2/BICseq_Norm/NBICseq-norm.pl
tmp=$norm/tmp
mkdir -p $tmp
readlen=100
stats=/home/minerva/WGS/fastq/$project/$patient/$library/BICseq2/$library.filt.stats.txt
fragsize=$(awk '/'average:'/ {print int($5);}' $stats)
out=$norm/$library.norm.out.$binsize.txt


#Options:
        #--help
        #-l=<int>: read length
        #-s=<int>: fragment size
        #-p=<float>: a subsample percentage: default 0.0002.
        #-b=<int>: bin the expected and observed as <int> bp bins; Default 100.
        #--gc_bin: if specified, report the GC-content in the bins
        #--NoMapBin: if specified, do NOT bin the reads according to the mappability
        #--bin_only: only bin the reads without normalization
        #--fig=<string>: plot the read count VS GC figure in the specified file (in pdf format)
        #--title=<string>: title of the figure
        #--tmp=<string>: the tmp directory;


#Generate normalized bin files
perl $BICSEQ_NORM -b=$bin -l=$readlen -s=$fragsize --gc_bin --NoMapBin --fig=$library --title=$library --tmp=$tmp $config $out


rm -rf $tmp


done < <(tail -n +2 /home/minerva/WGS/fastq/$project/$input)

