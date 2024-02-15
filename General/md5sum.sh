#!/usr/bin/env bash

# Sort files
find . -name \*fastq.gz | sort
find . -name \*fastq.gz | LC_ALL=C sort

# Check File Identity
find . -type f -exec md5sum {} +
find . -type f -exec md5sum {} + | LC_ALL=C sort #individual checksum
find . -type f -exec md5sum {} + | LC_ALL=C sort | md5sum #whoel directory checksum

md5sum * > checklist.chk        # Doesn't go down sub directories and generate checklist
md5sum -c checklist.chk   # runs through the list to check them! Use this!

# Check File Attributes
summary (){echo "$(stat -c '%y' "$1") $(md5sum "$1")"}
export -f summary
find . -type f -exec bash -c 'summary "$0"' {} \; | LC_ALL=C sort | md5sum

# Write to txt
find . -name \*fastq.gz -exec md5sum {} + | LC_ALL=C sort > ./nymd5sum.txt #individual *gz files
find . -type f -exec md5sum {} + | LC_ALL=C sort | md5sum > ./nychecklist.chk #whole directory

# Check Two md5sum files
awk 'FNR==NR{a[$1]=$1;next}{print $0,a[$1]?a[$2]:"NA"}' md5sum.txt ctg-md5.fastq.txt | grep NA | awk '{print $1,$2}' > ./md5sumdiff.txt


