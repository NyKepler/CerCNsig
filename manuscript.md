# CerCNsig
Generation of Copy Number Signatures for Cervical Samples from High Grade Serous Carcinoma Patients

## 1. Introduction

Ovarian cancer is a heterogeneous disease consisting of different histological subtypes with potentially different origins. High-grade serous ovarian carcinoma (HGSC) is the most common and aggressive form
of epithelial ovarian cancer, and prior studies have shown the p53 signature lesions and serous tubal intraepithelial carcinomas (STICs) in the fallopian tube are likely to be the common biologic origin and
precursors of HGSC and could be detected in the cervical specimens collected from 20 months to 6 years before the diagnosis 1,2. The ultimate goal of this study is to detect genome-wide copy number aberrations (CNAs) in asymptomatic cervical cytology specimens and generate unique tumorigenic signatures for early detection of HGSC.

## 2. Cohort Design
The first cohort contains HGSC, BRCA1/2 and benign patients who were operated in southern Sweden hospitals between 2015-2017. Their Blood, plasma, endometrial biopsy, tumor or tissue as well as a liquid-based cytology (LBC) sample were collected at the surgery and stored at -80°. Their archival cytology slides including PapSmear, liquid-based Pap or ThinPrep® were retrieved from pathology departments. The second cohort contains archival methanol-based ThinPrep® LBC, slides and tumor tissue from HGSC patients diagnosed between 2018-2020.

## 3. Methods and Material
DNA was extracted by Qiagen silica column or Maxwell magnetic bead-based method based on the sample type and performance. DNA was measured by Qubit fluorometric assays and quality-checked by Illumina TruSeq FFPE Library Prep QC Kit (qPCR). 1-10ng DNA was used to prepare next-generation-sequencing (NGS) library by QIAseq FX DNA library Prep kit. The shallow or low-pass whole genome sequencing (WGS) is performed on Illumina NovaSeq 6000 platforms.


#### Bioinformatic workflow
[BINP52_CNA_Framework](https://github.com/IngridHLab/BINP52_CNA_Framework), a pipeline to generate copy number profiles and detect copy number signatures from shallow whole genome sequening (sWGS) samples.

The version of tools and packages to be used will be specified in each step (see chapter 3). The scripts within the pipeline are based on Python (v3.11.6) and R (v4.3.2).
- (1) Preprocessing. This step includes quality assessment and quality trimming on the raw reads. (`Fastp` will be used for QC and trimming, together with `fastqc` and `multiQC` to generate the QC reports.)
- (2) Alignment. The human reference genome will be indexed. And the reads will be mapped to the reference genome. (`BWA` will be used for both indexing and alignment.)
- (3) Clean-up. After alignment, the SAM files will be sorted and the PCR duplicates will be marked and removed. Also, the .sorted.deduplicated.sam will be converted to BAM files. The BAM files will be indexed for later analysis. (`Picard` will be used for sorting SAM, marking duplicates, removing duplicates and converting SAM to BAM. `samtools` will be used for generating the clean_up stats and for indexing the BAM files.)
- (4) Relative copy number profile. The BAM files will be analyzed through fixed-size binning, filtering, correction, normalization to generate the read counts per bin. This data will then used for segmentation of bins and generating the relative copy number profile. (`QDNAseq` will be used for this step.)
- (5) Ploidy and cellularity solutions. The output file from `QDNAseq` contains relative copy number, and we need to estimate ploidy and cellularity in our samples to generate our final absolute copy number profile for comparison. (`Rascal` will be used for this step to find the solutions that best fit our study samples.)
- (6) Absolute copy number profile. We will further use other information (such as TP53 allele frequency) inferring the tumour fraction to select the best ploidy and cellularity solution. We apply this best solution to our relative copy number profile, and generate the final absolute copy number profile for each sample. (`Rascal` will be used for this step.)
- (7) Comparison with the recent HGSC signatures (n=7). The functions should be loaded from the github repository: https://bitbucket.org/britroc/cnsignatures.git .
- (8) Comparison with the Pan-Cancer signatures (n=17). The package `CINSignatureQuantification` will be used to generate the samply-by-component matrix for the Pan-Cancer chromosomal instability signatures.
- (9) Comparison with the panConusig signatures (n=25). Tools including `Battenberg` (`alleleCounter`, `impute2` and `beagle5` were included in this package), `ASCAT.sc` and `panConusig` will be used in this step.  
- (10) Generate and validate CerCN signatures (n=6) and the HGSC prediction model. Tools including mixture modeling `flexmix`, unsupervised clustering `NMF`, cosine similarity `lsa`, k-fold particitioning `caret` and `randomforest` will be used in this step.

![Figure_2_240705](https://github.com/NyKepler/CerCNsig/assets/111468388/28e80dea-300a-476a-bb3b-7b96a26514aa)

