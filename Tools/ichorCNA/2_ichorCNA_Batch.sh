#!/usr/bin/env bash

#conda activate R
#Execue script ./WGS/ichorCNA/Script/2_ichorCNA_Batch.sh


seq="CTG_2023_209"
project="MaNiLa_CS3"
input=MaNiLa_CS3_input.csv
group=MaNiLa_CS3_VS

# Create a txt file for ichorCNA parameters of all samples
version=v3
result=/home/minerva/WGS/ichorCNA/Data/MaNiLa/Parameters.$group.$version.txt
printf "Library\tTumorFrac\tPloidy\tSubclone\n" > $result

while IFS=, read -r patient type sample library fastq bin
do

pwd=/home/minerva/WGS/Data/$project/$patient
mkdir -p $pwd

#Create batch job for different binsizes
assay=$pwd/ichorCNA/$sample
mkdir -p $assay
cd $assay

wig=$library.$bin"kb".wig

pkg=/home/minerva/miniconda3/envs/R/share/r-ichorcna-0.3.2-1

# Remove old analysis plots to free space
rm -rf $assay/$library

# Benchmarking of ichorCNA using metastatic breast/prostate patient cfDNA and healthy donor cfDNA reveals a lower limit of sensitivity for detecting the presence of tumor to be 0.03 TFx (3%). That is, an ichorCNA-estimate of > 0.03 TFx will reliable indicate the presence of tumor for a cfDNA sample sequenced to ~0.1x whole genome coverage. An estimate of < 0.03 TFx indicates lowly detectable (below 3% but the estimate is less accurate) or absence of tumor-derived DNA. In the benchmark, at a 0.03 TFx cut-off, ichorCNA had a 91% specificity (correctly classified 20/22 healthy donors) and 95% sensitivity (classified 1125/1288 cancer patient mixtures). When the data quality differs, such as sequencing coverage or overall data variance, and a cancer type that is distinctly different than breast and prostate, manual inspection of cases near the 0.03 TFx cutoff and/or tuning of parameters is recommended.

# When selecting cfDNA samples for standard depths of whole exome sequencing (e.g. mean target coverage ~150x), we recommend samples with ichorCNA estimates of > 0.1 TFx (10%). From  the benchmarking, ichorCNA correctly estimated samples with > 0.1 TFx for 91% (606/669 patient-healthy donor mixtures) and correctly estimated samples with < 0.1 TFx for 96% (613/641 mixtures).

# see more: https://github.com/broadinstitute/ichorCNA/wiki/Interpreting-ichorCNA-results

# MAD value in the params.txt file? If it's > 0.15, then it may be too noisy

# bin size should not be smaller than 500kb if there is no matched normal sample. Ideal 1000kb is best for stand-alone samples.

# version 1: --id $library --WIG $assay/$wig --ploidy "c(2,3,4,5)" --normal "c(0, 0.05, 0.1, 0.2)" --maxCN 5 , --centromere $pkg/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt, --includeHOMD False, --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" --estimateNormal True --estimatePloidy True --estimateScPrevalence True --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --normalPanel $pkg/extdata/MaNiLa_Benign_PoN_$bin"kb"_median.rds --outDir $assay

# Vesion 2: --id $library --WIG $assay/$wig --NORMWIG $nassay/$nwig --ploidy "c(2,3,4,5)" --normal "c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)" --gcWig $pkg/extdata/gc_hg19_"$bin"kb.wig --mapWig $pkg/extdata/map_hg19_"$bin"kb.wig --centromere $pkg/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt --normalPanel $pkg/extdata/MaNiLa_Benign_"$type"_PoN_"$bin"kb_median.rds --outDir $assay


# Version 3: based on ACE and Rascal TF to fine tune the parameters
# Low tumor content samples (early stage disease): https://github.com/broadinstitute/ichorCNA/wiki/Parameter-tuning-and-settings
if [[ $type == *"TPVSs"* ]]; then 
Rscript $pkg/scripts/runIchorCNA.R --id $library.$bin"kb" --WIG $assay/$wig --ploidy "c(2,3)" --normal "c(0.85, 0.9, 0.95, 1)" --gcWig $pkg/extdata/gc_hg19_"$bin"kb.wig --mapWig $pkg/extdata/map_hg19_"$bin"kb.wig --centromere $pkg/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt --normalPanel $pkg/extdata/PoN/MaNiLa_Benign_ArchivalVS_PoN_"$bin"kb_median.rds --outDir $assay --plotFileType "png" --chrs "c(1:22)" --chrTrain "c(1:22)" --estimateScPrevalence False --scStates "c()" --maxCN 4

elif [[ $type == *"MaNiLaVS"* ]]; then 
Rscript $pkg/scripts/runIchorCNA.R --id $library.$bin"kb" --WIG $assay/$wig --ploidy "c(2,3)" --normal "c(0.9, 0.95, 1)" --gcWig $pkg/extdata/gc_hg19_"$bin"kb.wig --mapWig $pkg/extdata/map_hg19_"$bin"kb.wig --centromere $pkg/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt --normalPanel $pkg/extdata/PoN/MaNiLa_Benign_MaNiLaVS_PoN_"$bin"kb_median.rds --outDir $assay --plotFileType "png" --chrs "c(1:22)" --chrTrain "c(1:22)" --estimateScPrevalence False --scStates "c()" --maxCN 4

fi

#extract parameters from "params" files from all samples to one file
params=$assay/$library.$bin"kb".params.txt 
subclone=$(awk '/'Subclone'/ {print ($3);}' $params)
mad=$(awk '/'GC-Map'/ {print ($4);}' $params)
sed -n '2{p;q}' $params > tmp.txt
cat tmp.txt | sed s/$/"\t"$subclone/g | sed s/$/"\t"$mad/g >> $result
rm -f tmp.txt


done < <(tail -n +2 /home/minerva/WGS/ichorCNA/Data/MaNiLa/$input)










