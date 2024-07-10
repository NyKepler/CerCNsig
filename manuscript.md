# Cervical copy number signatures for early detection of high-grade serous tubo-ovarian carcinoma

## 1. Introduction
Ovarian cancer is a heterogeneous disease consisting of different histological subtypes with potentially different origins. High-grade serous ovarian carcinoma (HGSC) is the most common and aggressive form
of epithelial ovarian cancer, and prior studies have shown the p53 signature lesions and serous tubal intraepithelial carcinomas (STICs) in the fallopian tube are likely to be the common biologic origin and
precursors of HGSC and could be detected in the cervical specimens collected from 20 months to 6 years before the diagnosis 1,2. The ultimate goal of this study is to detect genome-wide copy number aberrations (CNAs) in asymptomatic cervical cytology specimens and generate unique tumorigenic signatures for early detection of HGSC.

## 2. Cohort Design
The first cohort contains HGSC, BRCA1/2 and benign patients who were operated in southern Sweden hospitals between 2015-2017. Their Blood, plasma, endometrial biopsy, tumor or tissue as well as a liquid-based cytology (LBC) sample were collected at the surgery and stored at -80°. Their archival cytology slides including PapSmear, liquid-based Pap or ThinPrep® were retrieved from pathology departments. The second cohort contains archival methanol-based ThinPrep® LBC, slides and tumor tissue from HGSC patients diagnosed between 2018-2020.

## 3. Copy number analysis
DNA was extracted by Qiagen silica column or Maxwell magnetic bead-based method based on the sample type and performance. DNA was measured by Qubit fluorometric assays and quality-checked by Illumina TruSeq FFPE Library Prep QC Kit (qPCR). 1-10ng DNA was used to prepare next-generation-sequencing (NGS) library by QIAseq FX DNA library Prep kit. The shallow or low-pass whole genome sequencing (WGS) is performed on Illumina NovaSeq 6000 platforms.


### Bioinformatic workflow
[BINP52_CNA_Framework](https://github.com/IngridHLab/BINP52_CNA_Framework), a pipeline to generate copy number profiles and detect copy number signatures from shallow whole genome sequening (sWGS) samples.

[CerCNsig](https://github.com/IngridHLab/CerCNsig), a workflow to optimal absolute copy number profiles and generate of copy number signatures for cervical Samples from high grade serous carcinoma patients.

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

### Operation guide

#### 1. Workflow *Part I Solutions* in [BINP52_CNA_Framework](https://github.com/IngridHLab/BINP52_CNA_Framework)
```
# activate the conda environment installing snakemake
conda activate snakemake
# move into the working directory where you downloaded the github repository
cd <path to the downloaded github repository>
# run the pipeline, for example, 30 threads
snakemake --use-conda --configfile config/config.yaml --cores 30 --snakefile workflow/Snakefile_solution.smk
```
#### 2. Workflow *Part II Signatures* in [BINP52_CNA_Framework](https://github.com/IngridHLab/BINP52_CNA_Framework)
Note: this step works for HGSC CN signatures as well as pan-cancer CIN signatures, but not includes the panConusig signatures.
```
# activate the conda environment installing snakemake
conda activate snakemake
# move into the working directory where you downloaded the github repository
cd <path to the downloaded github repository>
# run the pipeline, at least 30 threads is suggested
snakemake --use-conda --configfile config/config.yaml --cores 30 --snakefile workflow/Snakefile_CNsig.smk
```
#### 3. Workflow *Part III panConusig* in [BINP52_CNA_Framework](https://github.com/IngridHLab/BINP52_CNA_Framework)
Note: for panConusig, it requires "normal-tumor" pair, the sample sheet is different from other signatures (see example `config/sample_pair_3.tsv`). Furthermore, it starts again from the BAM files to generate allele-specific copy number profiles which are completely different from the other two signatures, so it is good to run separately. Besides, it would be better to run in local conda environment instead of within the Snakemake pipeline to aviod dependency conflicts (such as Java) as well as pathway errors of reference files.  
Please see suggestion on setting up the environment in [section 3.3](https://github.com/IngridHLab/BINP52_CNA_Framework/blob/main/README.md).
```
# Firstly, we need to add our local true working directory in the impute reference files
cat resources/battenberg/battenberg_impute_v3/impute_info00.txt | sed 's#<path_to_impute_reference_files>#resources/battenberg/battenberg_impute_v3#g' > resources/battenberg/battenberg_impute_v3/impute_info.txt

# Secondly, we need to build up the environment with required dependencies
Rscript workflow/scripts/panConusig_env.R

# Then, run panConusig generation step by step.
# Please remember to change the directories in the scripts to adapt to the user situation
Rscript workflow/scripts/panConusig_pair_local_1.R
Rscript workflow/scripts/panConusig_pair_local_2.R
Rscript workflow/scripts/panConusig_pair_local_3.R
```
#### 4. Workflow *Part VI [CerCNsig](https://github.com/IngridHLab/CerCNsig/tree/main/Tools/CNsignatures)* 
The [CerCNsig_filt](https://github.com/NyKepler/CerCNsig/tree/main/Tools/CNsignatures) version 2024 based on the absolution copy number profiles of HGSC Cervical samples generated from the https://github.com/IngridHLab/BINP52_CNA_Framework pipeline. Cervical samples were selected based on their HGSC CN signatures in Macintyre et al. 2018 https://github.com/markowetzlab/CNsignatures: samples with similiarity more than the first three signatures (S1-S3). Those cervical samples were considered to be CNA enriched instead of filtering the cervical samples using the cellularity from ACE/Rascal/ichorCNA estimation and mauanlly inspection which could be not completely accurate.

```
#' Select the cervical samples by removing the samples with noisy copy number profiles based on the QDNAseq observed std and estimated have HGSC CNsig > 3 (42 samples remained).

#' Compare underlying distribution of six HGSC CNsignature copy number features in tumor and all cervical samples.
Rscript Tools/CNsignatures/Plot_Copy_Number_Features_Distribution.R

#' Mixture modelling and extract 36 unique components from the selected 42 HGSC cervical samples.
./Tools/CNsignatures/Extract_VS_Sample_by_Component_filt.sh
Rscript Tools/CNsignatures/Extract_VS_Sample_by_Component_filt.R

#' Optional: mixture modelling to identify unique components otherwise the script above is used to extract sample_by_component matrix based on the predifined components:
type <- "HGSC_VS_good_CNsig_over_3"
CN_exposure_over_3 <- readRDS("/home/researcher/CerCNsig/Randomization/CN_exposure_over_3.rds")
input <-inner_join(input, CN_exposure_over_3)
write.table(input, "Sample.List.rascal.HGSC_VS_good.filt.txt", sep ="\t", quote = F, row.names = F, col.names = T)
VS_components <- fitMixtureModels(CN_features)
saveRDS(VS_components, file = paste0(type, ".Components.rds"))

#' Generate new version of CerCN signatures (6 sigs) from the selected 42 HGSC cervical samples.
Rscript Tools/CNsignatures/CerCNsig.R

#' Compare the cosine similarity of each CerCN signatures across Benign, RRSO and HGSC groups.
#' Compare the cosine similarity of CerCN signatures of different diagnostic time points within the same group.
Rscript Tools/CNsignatures/CerCNsig_Cosine_Matrix.R
Rscript Tools/CNsignatures/Plot_Signature_Cosine_Comparison.R

#' Defining the CerCN signatures by comparing the component weights of each signature to the HGSC CNsignatures: a histogram containing the relative weighting of each of the components, colour coded as the feature distribution.
Rscript Tools/CNsignatures/Plot_Signature_Component_weights.R
```
#### 5. Workflow *Part V [RandomForest](https://github.com/IngridHLab/CerCNsig/tree/main/Tools/RandomForest) classification model*
To apply random forest modeling to generate a prediction or classification model on cervical samples, we first filter away noisy cervical samples based on QDNAseq CN profile and the differences between expected standard deviation (Eσ) and measured standard deviation, and then applied all 36 components from 163 selected cervical samples from both HGSC and benign patients in the model.

```{bash randomforest k-fold validation}
#' Generate and validate randomforest model on merged sample_by_component (163 x 36) matrix of HGSC and benign.
Rscript Tools/RandomForest/Randomforest_CV_CerCNsig_filt.R
```
