#!/usr/bin/env bash

# This is the 4th step of the BICseq workflow
# conda activate R
# execue by ./WGS/BICseq2/Script/4_BICseq_Segment.sh

project="MaNiLa_CS3"
input=MaNiLa_CS3_input.csv

while IFS=, read -r patient type sample library fastq bin
do

pwd=/home/minerva/WGS/Data/$project/$patient/BICseq2/$sample
mkdir -p $pwd


bin=10 #Default 100bp but according to the 2016 publication We normalized data and set the initial bin size as 10 bp for segmentation. The penalty parameter λ was chosen as 1.2
binsize="$bin"bp
lambda="1.2" #The smoothness of the CNV profile. The larger the parameter is, the less segments the file profile will have. Default is 2. Based on the latest publication the penalty parameter λ was chosen as 1.2. 1.5 was used in CS1
seg=$pwd/Seg_Out_"$binsize"_"$lambda"
mkdir -p $seg

norm=$pwd/Norm_Out_$binsize

#Create config file with header and fill in location of each file into config file
config=$seg/$library.seg-config.$binsize.txt
printf "chromName\tbinFileNorm\n" > $config #Create config file with header for stand-alone sample!

for x in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22} 
do 
chromName=chr"$x"
binFileNorm=$norm/$library.chr$x.$binsize.norm.bin 
printf "$chromName\t$binFileNorm\n" >> $config #Stand-alone control samples!!
echo $x
done

#Segmentation on sample without control (e.g.,blood or normal tissue) using BICseq2-seg: 
BICSEQ_SEG=/home/minerva/WGS/WGS_v2/BICseq2/BICseq_Seg/NBICseq-seg.pl
tmp=$seg/temp
mkdir -p $tmp
png=$seg/$library.CNV.$binsize.png
out=$seg/$library.CNV.$binsize.txt


	#--lambda=<float>: the (positive) penalty used for BICseq2
	#--tmp=<string>: the tmp directory
        #--help: print this message
        #--fig=<string>: plot the CNV profile in a png file
        #--title=<string>: the title of the figure
        #--nrm: do not remove likely germline CNVs (with a matched normal) or segments with bad mappability (without a matched normal)
        #--bootstrap: perform bootstrap test to assign confidence (only for one sample case)
        #--noscale: do not automatically adjust the lambda parameter according to the noise level in the data
        #--strict: if specified, use a more stringent method to ajust the lambda parameter
        #--control: the data has a control genome
        #--detail: if specified, print the detailed segmentation result (for multiSample only)



 
perl $BICSEQ_SEG --tmp=$tmp --lambda=$lambda --detail --bootstrap --fig=$png --title=$library.$binsize.$lambda $config $out #For stand-alone sample without control!!
rm -rf $tmp



done < <(tail -n +2 /home/minerva/WGS/fastq/$project/$input)



