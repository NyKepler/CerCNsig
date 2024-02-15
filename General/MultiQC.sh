#!/usr/bin/env bash

patientID="MaNiLa_HGSOC_137"

cd /home/minerva/WGS/Data/MaNiLa_CS1/"$patientID"_210911
for i in {6,7,8,9,29,35}
do
sampleID="CS1_"$i"_S"$i
cd $sampleID

multiqc .

cd ..

echo $i
done
