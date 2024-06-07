# CerCNsig
Generation of Copy Number Signatures for Cervical Samples from High Grade Serous Carcinoma Patients

## 1. Introduction

Ovarian cancer is a heterogeneous disease consisting of different
histological subtypes with potentially different origins. High-grade
serous ovarian carcinoma (HGSC) is the most common and aggressive form
of epithelial ovarian cancer, and prior studies have shown the p53
signature lesions and serous tubal intraepithelial carcinomas (STICs) in
the fallopian tube are likely to be the common biologic origin and
precursors of HGSC and could be detected in the cervical specimens
collected from 20 months to 6 years before the diagnosis 1,2. The
ultimate goal of this study is to detect genome-wide copy number
aberrations (CNAs) in asymptomatic cervical cytology specimens and
generate unique tumorigenic signatures for early detection of HGSC.

## 2. Cohort Design

The first cohort contains HGSC, BRCA1/2 and benign patients who were
operated in southern Sweden hospitals between 2015-2017. Their Blood,
plasma, endometrial biopsy, tumor or tissue as well as a liquid-based
cytology (LBC) sample were collected at the surgery and stored at -80°.
Their archival cytology slides including PapSmear, liquid-based Pap or
ThinPrep® were retrieved from pathology departments. The second cohort
contains archival methanol-based ThinPrep® LBC, slides and tumor tissue
from HGSC patients diagnosed between 2018-2020. Till now total 320 DNA
samples have been processed and analysed in which 7 patients are Benign,
4 BRCA and 19 HGSC upon that 7 HGSC are from cohort 2 and rest from
cohort 1.

## 3. Methods and Material

DNA was extracted by Qiagen silica column or Maxwell magnetic bead-based
method based on the sample type and performance. DNA was measured by
Qubit fluorometric assays and quality-checked by Illumina TruSeq FFPE
Library Prep QC Kit (qPCR). 1-10ng DNA was used to prepare
next-generation-sequencing (NGS) library by QIAseq FX DNA library Prep
kit. The shallow or low-pass whole genome sequencing (WGS) is performed
on Illumina NovaSeq 6000 platforms. The demultiplexed FastQ data will be
processed with BWA (alignment), SAMtools (cleanup), Picardtools
(duplicate) and gatk (BQSR). Copy number analysis is performed by
modified QDNAseq, BICseq2 and HMMcopy (prior ichorCNA). Tumor fraction
(TF) and ploidy is estimated by ACE, ichorCNA or Rascal. Copy number
signatures are obtained by CNsignatures and focal CNA by GISTIC2. Other
data analysis and visualization are done in R and RStudio.

<p align ="center">
<img width="500" alt="image" src="https://github.com/NyKepler/CerCNsig/assets/111468388/4ecd3f0f-4920-4130-aae5-8e6d4fc154b3">
</p>

## 4. Data Analysis and Results

This report is mainly focusing on WGS data analysis procedure and
further validation is required to fully understand the cancer biology
and mechanism behind the findings.

#### 1. Sequencing Data Overview

Fastq files were aligned by *BWA* with different parameters recommended
by the downstream program to generate alignment BAM files. The BAM files
then sorted and indexed by *Samtools* and alignment metrics were
extracted by *Samtools stats* and *Samtools Flagstat*. Duplicates and
other unwanted reads can be marked and removed by either *Samtools* or
*Picard*. Those results from all samples were aggregated into one html
report by *MultiQC* for easy access.

In the POC study we got 36-63 million read-pairs (mrp) with 20-34%
duplication and 28-47 mrp after removing duplicates. The coverage is
0.9-1.6X.

In the cohort study we got 23-101 mrp with 22-78% duplication and 13-69
mrp after deduplication. The coverage is 0.4-2.2X. The coverage is
calculated based on the formula below:

$$
C=\frac{L*N}{G}
$$

-	G is the haploid genome length, 3.2\^
-	L is the read length, ex. 100bp
-	N is the number of reads, 30M
-	C stands for coverage, be roughly 1X

#### 2. [QDNAseq](https://github.com/NyKepler/CerCNsig/tree/main/Tools/QDNAseq) Copy Number Analysis 

Description from the Bioconductor package "QDNAseq" developed by Ilari
Scheinin and others (1): Quantitative DNA sequencing for chromosomal
aberrations. The genome is divided into non-overlapping fixed-sized
bins, number of sequence reads in each counted, adjusted with a
simultaneous two-dimensional loess correction for sequence mappability
and GC content, and filtered to remove spurious regions in the genome.
Downstream steps of segmentation and calling are also implemented via
packages DNAcopy and CGHcall, respectively.

Following the instruction of the 'Introduction to QDNAseq' in the
package, in-house bin annotations for the hg19 build of human reference
genome were generated in order to optimize the performance of CN
analysis on 100bp paired-end sequencing setup which differs from the
original 50bp single end. A new control set with 58 human genomes is
obtained from the 1000 Genomes Project Phase 3 using the following
filters:

    Keep cases that are on Illumina platform, low coverage, not withdrawn and paired-end sequenced in a single lane:
    - g1k <- g1k[g1k$INSTRUMENT_PLATFORM == "ILLUMINA", ]
    - g1k <- g1k[g1k$INSTRUMENT_MODEL == "Illumina HiSeq 2000", ]
    - g1k <- g1k[g1k$ANALYSIS_GROUP == "low coverage", ]
    - g1k <- g1k[g1k$LIBRARY_LAYOUT == "PAIRED", ]
    - g1k <- g1k[g1k$WITHDRAWN == 0, ]
    - g1k <- g1k[!g1k$PAIRED_FASTQ == "", ]

The reference genome was downloaded from the 1000 Genomes project and
indexed by BWA as to be used in the BWA alignment of the control
dataset. Reference genome should be indexed by picard and samtools if
using the GATK tool kits for removing duplication of the actual sample
dataset. Mappability file in bigWig format and blacklist files in bed
format was downloaded from UCSC accordingly.

General statistics, percent mapped and alignment metrics from the result
of *Samtools stats* were summarized by MultiQC into one single html
report for easy comparison. One sample named SRR400036 was excluded due
to higher error rate and less reads than others.

Locate the bam files and run the Rscript below to create new annotations
on all 9 bin sizes (1, 5, 10, 15, 30, 50, 100, 500 and 1000kb). Note:
the median residuals of the LOESS fit of the new annotations don't
correlate well to the pre-generated annotations (under 0.5) due to
differences on sequencing read length (2x100 vs 50bp) and depth (\~2 vs
\~0.1X) of the two control datasets. Note, this should not be run in
markdown mode due to heavy computation.

Spot noisy samples: https://github.com/ccagc/QDNAseq/issues/35

#### 3. [ACE](https://github.com/NyKepler/CerCNsig/tree/main/Tools/ACE) absolute copy number analysis: 
-    Ploidy and cellularity estimation as well as absolute copy number fitting on QDNAseq segmented relative copy number data. 

Absolute copy number estimation (ACE) scales relative copy number
signals from chromosomal segments to optimally fit absolute copy
numbers, without the requirement of additional genetic information like
SNP data. ACE derives an estimate of tumor purity together with ploidy.
ACE will calculate how well the segments fit (the relative error) to
integer copy numbers for each percentage of "tumor cells" (cells with
divergent segments). Note that it does not estimate for a lower
percentage than 5 (0.05). The squaremodel function chooses the ploidy
and cellularity (also called tumor fraction or tumor purity) on the
median bin segment value (standard) with the lowest relative error,
which has been validated on a dataset of 253 ovarian carcinoma samples
(Brentons CNsignature study) and the ACE estimated tumor fraction from
copy number data showed higher correlation to the high-depth mutation
data (golden standard) than the model ABSOLUTE. Alternatively, the
runACE function allower user manually choose standard and ploidy to
obtain optimal cellularity.

The underlies mathematics and methodologies of this model are summarized
below. Segment data (signal intensity and number of bins) are obtained
from the QDNAseq-object using the R function rle(). For any
"cellularity" between 0.05 and 1 with increments of 0.01, the expected
relative signal is calculated for all integer copy numbers "copies"
between 1 and 12 as below:

$$
Signal^{expected} = Standard^{observed} * {\frac{cellularity * copies+2*(1-cellularity)}{cellularity * ploidy+2*(1-cellularity)}}
$$ 


The "standard" is the median segment value of all bins as mentioned before and ploidy is the general ploidy of the aberrant cells.

The difference between the segment value and the closest expected signal
is the "error" of the segment. The default setting calculates the error
of the fit as the root mean square error (RMSE) of all segments. To
account for segment length, segment errors are repeated as many times as
the number of bins the segment comprises. Alternatively, when mean
absolute error (MAE) is chosen as error method, the error of the fit is
the mean of the absolute segment errors (again segment length is
accounted for as described above).

$$
RMSD = \sqrt{\frac{\sum_{n=1}^{N}{(Signal^{expected}-Standard^{observed})}^2}{N}}
$$ 

$$
MAE = \frac{\sum_{n=1}^{N}{|Signal^{expected}-Standard^{observed}|}}{N}
$$ 

Additional, optional parameters like "penalty" (0 to 1) and
"penploidy" (0 to 1) to penalize for fits at low cellularity or ploidies
different from 2N can be specified by user. Default settings are 0.5 for
both "penalty" and "penploidy" but I have run some comparison to find
the better settings for our samples (see below).

$$
error = \frac{error}{cellularity^{penalty}}
$$

$$
error = error * (1 + |ploidy – 2|)^{penploidy}
$$ 

Screening for the parameters like penalty, penploidy in the model as
well as amount of read and bin sizes (30, 50 and 100kb) on Benign
patient Cervical samples to determine false positive limit and on tumor
samples for false negative cellularity threshold. According to the ACE
article, false positive results are highly unlikely using a penalty
factor of 0.1 on high quality, fresh, healthy control DNA sample.
However, penalty factor for optimal sensitivity should be trained using
the appropriate control in this case cytologic samples.

In order to optimize/reduce the penalty factor and increase the
sensitivity of tumor fraction estimation without causing false positive,
I run ACE squaremodel on QDNAseq segment results from 3 different
binsizes (30, 50 and 100 kb) with penalties from 0 to 1 (no to max
penalty) with 0.01 interval and penploidy always at 0.5 on benign
patient samples. The cellularities and ploidies were tested in two
different settings: version 1 is adapted to ovarian cancer and version 2
is default setting in the squaremodel. The model returns combinations of
cellularity and ploidy with minimum error as minima in the "minimadf"
list and order by error. The top minima with lowest error from each
squaremodel run is considered as the "bestfit". Note that an estimate of
1 (100% tumor cell percentage) by ACE generally indicates a negative
result or 'flat' absolute copy number profiles without CNAs. Any minima
with cellularity equal to 1 was considered as 'negative call' and saved
for the analysis below.

The optimization has been run on 6 sample types blood, plasma, tissue,
endome, MaNiLa Cervical and archival Cervical which indicate different coverage
depths and , 3 bin sizes during QDNAseq read counting and 2 distance
error methods in ACE. This have been run in two slightly different
versions as mentioned in the code above: version 1 is has been used for
validation of ACE on the Brentons CNsignature ovarian cancer cohort with
cellularity 1 - 100; ploidy 1.8 - 4.3, while version 2 is the default
setting with cellularity 5 - 100; ploidy 1 - 5. Note: Benign patient_124 
ovarian tissue (CS2_49) shows tumor content in all fittings due to
a distinct copy number gain of the entire chromosome 12, so her tissue
sample is excluded from the rest of analysis.

Each data point in the NCplots shows the penalty value with cellularity
= 1 (no CNA or no false positive). The line connects the estimating
penalty value at the upper bound of 95% confident interval for different
sample types from benign patients, and color of the line indicates the
error method during modeling. According to benchmark result from ACE
paper, MAE error method tends to produce slightly more accurate results,
and penalty can increase accuracy but reduce sensitivity when
cellularity or tumor fraction is low. I ran the lowest penalties as well
as the more sensitive binsize with lowest penalty (depending on the
genomic size and abundance of CNAs a larger bin size improves
sensitivity but may decrease accuracy) to benchmark the performance on
MaNiLa samples.

The top fitting with smallest error is chosen for all samples in the
first place but changes were made for individuals after manually
inspection of matrixplots. By default, ACE does not estimate lower
percentages than 5%.

Postanalysis based on the chosen fitting (cellularity and ploidy) with
given standard = 1 is performed to scale bin values and segment values
to absolute copies using the formula below, and a segment table with
absolute copy number is generated using postanalysisloop function: 
	
$$
copies = ploidy + (signal − standard) \frac{cellularity * ploidy + 2 (1 − cellularity)}{cellularity * standard}
$$

To summarize the copy number alteration in our HGSC sample, over 15% of
them have amplification in regions containing MECOM, PIK3CA, KRAS, CCNE1
and AKT2, while 15% of them have deletion in CDKN2A/B, NF1, ERBB2
(proximate) and CCNE1. Large segment amplification on 3q, 10p, 12p and
20q, and deletion on 4p, 8p, 11p and 16q. According to the Z.Cheng's
paper from McNeish lab, amplification in genes like PTEN, CCNE1, CCND1,
AKT2 and deletion in TP53 are more frequent in early stage HGSC, while
amplifiation in PTEN and KRAS are more seen in later stage. More detail
focal CNA analysis will be done by GISTIC2.

#### 4. [Rascal](https://github.com/NyKepler/CerCNsig/tree/main/Tools/Rascal) absolute copy number analysis: 
-    Benchmark with ACE shows super for HGSC tumor samples but not ideal for cytology Cervical. 

*Rascal* was developed based on concepts of ACE and further adapted to
the biological nature of ovarian tumor samples. The mathematics of this
approach assume a single dominant clone and not suitable for very
heterogeneous tumor samples. The two distance functions MAD (the mean
absolute deviation) and RMSD (the root mean square deviation) are based
on the difference between expected copy numbers are to their closest
whole integers, which unlike ACE squaremodel the MAE or RMSE is based on
difference between expected copy numbers and the standard 1. According
the validation in the Rascal paper, MAD has better performance (higher
accuracy) in recapitulating the cellularity and ploidy than RMSD in
ovarian cancer cell lines and PDX model. Moreover, Rascal doesn't
utilize penalty factors like ACE to penalize fittings with low
cellularity or ploidies that diverge from the diploid 2N state, and
penalty on ploidy is unsuitable for cancer types where whole genome
duplication is frequent. Note, I didn't optimize the penploidy in ACE
and just keep it moderate at 0.5 as the original ACE validation on
Brenton's OC cohort.

The model identifies best fit solutions based on the lowest MAD distance
value, and while in cases where multiple competing solutions with the
lowest MAD distance are available, solution with the higher cellularity
is prefered if tumor purity is known.

For clinical sample like tumor tissue, the uncertainty of tumor purity
makes it difficult to choose absolute copy number fitting and
distinguish between competing best fits. Additional information like SNV
allele fractions or an accurate estimate of the tumor cellularity from a
pathologist, would be required to obtain reliable ACN fits. For example,
comparing the expected TP53 mutant allele fractions (MAFs) to the
empirical TP53 MAF measured from Tagged-Amplicon sequencing (TAMseq)
before entering an annotation and solution refinement stage. The
expected TP53 MAF is calculated using the formulation below:

$$
expectedTP53MAF= \frac{copies_{TP53} * cellularity}{copies_{TP53} * cellularity + 2(1-cellularity)}
$$

The default setting in the'fit_absolute_copy_number' function is listed
below and unlike the original ACE penalty or penploidy. My guess is that
they have validated those parameters on their OC cell lines, PDX models
and large HGSC patient cohort. In order to make a fair comparison, I
process the HGSC tumor segment data in 100kb from QDNAseq with some
modifications on max/min ploidy and min cellularity to keep it similar
to the ACE squaremodel version 1. The distance method is kept as default
MAD as well. Additional model adjusting was done on pat_48, pat_81 and
pat_118 using their TP53 MAF measured by NGS which also validated by
IBSAFE in the previous MaNiLa study. The precise TP53 copy number status
can only be obtained from the QDNAseq segmentation result at 30kb but
not at 100kb, however one can use the value from a much wider segment if
100kb is prefered. The analysis I demonstrate below is performed on the
30kb segmented data.

    --min-ploidy=1.8 The minimum ploidy to consider (default: 1.25)

    --max-ploidy=4.3
        The maximum ploidy to consider (default: 5.25)

    --min-cellularity=0.01
        The minimum cellularity to consider (default: 0.2)

    --max-cellularity=MAX-CELLULARITY
        The maximum cellularity to consider (default: 1)

The TP53 expected MAF can be calculated and used to adjust the solutions
previously, and the new solution which gives a 'TP53 tumor fraction'
closest to the observed allele fraction which is NGS_MAF or IBSAFE_MAF
in our study, would be most accurate solution. According to Rascal
paper, when the TP53 MAF difference (MAFdiff) between the expected and
the observed is bigger than 0.3 it means rascal failed to find a matched
solution, and if observed TP53 MAF is not available or below 0.3 (low
cellularity/purity) the fitting solutions should be manually inspected.

#### 5. [ichorCNA](https://github.com/NyKepler/CerCNsig/tree/main/Tools/ichorCNA) copy number analysis: 
-    What is the limit detect tumor fraction in WGS data. 

*ichorCNA* was originally designed for estimating circulating tumor DNA
in the plasma using ultra-low-pass whole genome sequencing (0.1x
coverage). In contrast to the previous two methods which perform
analysis on the output of the QDNAseq, the ichorCNA algorithm uses a
hidden Markov model (HMM) for the probablisitc modeling and works in
sequencing coverage down to 0.1X. The read count file in WIG format is
created from the indexed BAM file with a desired genomic window (ex.
50kb or 1000kb) by HMMcopy suite readCounter.

Alternative, one can convert the normalized read counts from QDNAseq to
WIG format and perform the segmentation in ichorCNA by HMMcopy. The wig
files I used for the analysis however were generated by HMMcopy
readCounter and should analyse the QDNAseq read counts to see if any
better.

*ichorCNA* also utilizes panel of normal or normal reference (ex. blood)
as well as pre-defined fraction of normal contamination during model
fitting to increase the specificity and sensitivity. Since we have
different sample types which were preserved either fresh frozen or fixed
and they were prepared respectfully for the whole genome sequencing, so
I generated panel of normal for each individual sample type using benign
patients blood, plasma, endome (endometrium biopsy), MaNiLaVS,
archivalVS, tissue and FFPE. Those PoNs should help to correct
systematic biases arising from library construction, sequencing platform
and DNA-specific artifacts. Note: 6 benign samples were excluded in the
PoNs as they have somewhere above 10 percent cellularity in the ACE
data, including pat_124 tissue which has aberration on the entire
chromosome 12. I have to address that the pipeline seems not using PoN
when there is normal reference available.

    #' Before start: Create PoN in R
    outfile=MaNiLa_Benign_"$type"_PoN_"$bin"kb
    pkg=/home/minerva/miniconda3/envs/R/share/r-ichorcna-0.3.2-1

    Rscript $pkg/scripts/createPanelOfNormals.R 
      --filelist $pwd/$wiglist 
      --gcWig $pkg/extdata/gc_hg19_"$bin"kb.wig 
      --mapWig $pkg/extdata/map_hg19_"$bin"kb.wig 
      --centromere $pkg/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt 
      --outfile $outfile

Based on the ploidy and cellularity range estimated by ACE together with
Rascal, I modified the default settings (according to
<https://github.com/broadinstitute/ichorCNA/wiki/Parameter-tuning-and-settings>)
to improve the estimation especially for those have cellularity below
0.5 and it is recommended to have at least 1X coverage for those
samples. According to the previous summary on QDNAseq used readcount,
most of our samples have over 1X except blood and plasma have about 0.5X
coverage.

    1. Initialize tumor fraction parameter: 
    --normal "c(0.95, 0.99, 0.995, 0.999)"
    Initialize the non-tumor (1 minus tumor fraction) to expected values, such as 5%, 1%, 0.5%, 0.1%. ichorCNA will still estimate the tumor fraction but having these initial starting values can help the EM step find better global optima.

    2. Set initial ploidy to diploid
    --ploidy "c(2)"
    It will be difficult to predict the ploidy value for low tumor fraction cases.

    3. Reduce number of copy number states
    --maxCN 4
    Reducing the state space will help reduce complexity. If you know from a prior sample (e.g. tumor biopsy) that there are large high level copy number events, you can set this to 4.

    4. Do not account for subclonal copy number events
    --estimateScPrevalence FALSE --scStates "c()"
    Subclonal events are difficult to detect for low tumor fraction, these can we turned off.

    5. Train and analyze autosomes only
    --chrs "c(1:22)" --chrTrain "c(1:22)"
    Exclude chrX in the analysis and training to reduce complexity.

Next we performed the ichorCNA version 3 fitting with additional fine
tuned parameters based on ACE and Rascal tumor fraction results, and the
results will be compared and used to further adjust the ACE fittings.
The final fittings will be used to extract the absolute segmented copy
number.

It is recommended to check the GC-Map correction MAD value to determine
the data quality. A MAD value greater than 0.3 usually indicates a high
variance in the data and it is likely too noisy. Normally GC and
mappability biases are corrected during analysis but a high MAD may
indicate noise or artifact that is not being corrected. On the other
hand, a MAD value less than 0.15 with coverage at least 0.1X are
considered to be high quality data.

The benchmarking of ichorCNA with only 0.1X coverage, using metastatic
breast/prostate patient cfDNA and healthy donor cfDNA, shows cellularity
\> 0.03 will reliable indicate the presence of tumor for a cfDNA sample
, while \< 0.03 indicates low (less reliable) or absence of
tumor-derived DNA. With 0.03 cut-off, ichorCNA had a 91% specificity and
95% sensitivity, but when the data quality differs or a distinctly
different cancer type, manual inspection of cases or additional tuning
of parameters is recommended.

The final estimation is revised manually that the tumor samples results are based on Squaremodel or Rascal while the rest are mostly based on ichorCNA with 0.03 as threshold meaning anything below that will be called 0 cellularity and normal diploidy.

#### 6. Copy Number Signatures Analysis: 
-    Whole genome copy number signature analysis to extract HGSC related signatures from pre-diagnotic/pre-sytompatic cytology samples. 

According to *CNsigntures* the analysis just require a list of absolute
segmented copy number tables which contains "chromosome", "start",
"end", "segVal".. That seems easy but no that simple as you think! In
the original validation for HGSC, samples were first analysed with sWGS
and TP53 Tagged amplicon seq, and then processed using ABCEL (might be
released as Rascal now) to generate absolute copy-number profiles. All
samples were manually inspected and the copy-number fit adjusted if
necessary (based on evidence from other samples from the same patient).
Following absolute copy-number fitting, the samples were rated using a
1-3 star system.

Quote: \* 1 star samples showed a noisy copy-number profile and were
considered likely to have incorrect segments and missing calls. These
were excluded from further analysis. \*\* 2 star samples showed a
reasonable copy-number profile with only a small number of miscalled
segments. These samples were used (with caution) for some subsequent
analyses. \*\*\* 3 star samples showed a high-quality copy-number
profile that was used in all downstream analyses.

Maybe I have spent too much an effort on the ACN estimation by trying
different tools and optimizing their performance, but combining the
results from ACE, Rascal and ichorCNA we should have a better estimation
of ACN profile on our low purity MaNiLa Cervical and archival Cervical samples and
be able to select the ones with true signal.

To better understand the principle and workflow of the *CNsignature* and
potentially generate our own sigatures for cytology samples, I will try
to demonstrate each step of the analysis with the codes in
'main_functions.R' and source the rest from the 'helper_functions.R'
instead just running the function:
        
    Before start: Prepare segment tables:

The ACE postanalysis provides a segment table with estimated ACN based
on the ploidy and cellularity fitting, and the column 'Chromosome',
'Start', 'End' and 'Segment_Mean2' are subtracted and renamed for the
analysis. Samples were filtered by cellularity not = 0 and number of
segments with mean value not = ploidy (indicating CNA) in order to
capture real signal in the signature analysis.

    Step 1. Extract 6 fundamental CN features:
    -   Segment size
    -   Breakpoint number (per ten megabases)
    -   Length of segments with oscillating copy-number
    -   Breakpoint number (per chromosome arm)
    -   Change-point copy-number
    -   Segment copy-number
    
    Step 2. (Optional) Apply mixture modelling to breakdown each copy number feature distribution into mixtures of Gaussian or mixtures of Poisson distributions using the [flexmix](https://cran.r-project.org/web/packages/flexmix/index.html) package. For HGSC analysis, the mixture model component definitions from the publication should be used and we skip this step. Otherwise, the optimal number of categories in each feature may differ between cancer types and should be defined by a mixed effects model.

    Step 3: Generated a sample-by-component matrix representing the sum of posterior probabilities of each copy-number event being assigned to each component. If the all_components argument is specified, then the sum-of-posteriors is calculated using these components, otherwise the component definitions from the publication (component_parameters.rds) are used.

Heatmap for matrix of 14 tumor samples by 36 components from the
validation versus matrix with 24 components found by *fitMixtureModels*.
The heatmap with 36 components are less colorful as the validation
cohort which contain 117 samples, which indicates our samples might be
more unified. When the components reduce to 24, the heatmap shows more
clusters between samples.

    Step 4: (Advanced) identify number of signatures using Non-negative matrix factorization [NMF package](https://cran.r-project.org/web/packages/NMF/index.html), which factorize the sample-by-component matrix into a signature-by-sample matrix and component by signature matrix. NMF demands the upper bound on the number of signatures be \<\< than both the sample and component number, so a signature search interval of [3,12] is used. The matrix factorization is run 1000 times for each number of signatures, each with a different random seed, and compared the cophentic, dispersion, silhouette, and sparseness scores for the signature-feature matrix (basis), patient-signature matrix (coefficients) and a consensus matrix of patient-by-patient across the 1000 runs. In addition we performed 1000 random shuffles of the input matrix to get a null estimate (Randomised) of each of the scores.

A value (number of signatures) defines the point of stability in the
cophenetic, dispersion and silhouette coefficients, and is the maximum
sparsity achievable above the null model (when observed still \>\> than
randomised) in the Basis (the red lines). According to that, the optimal
number of signatures would be 6. A signature number higher than this
would force the sparseness in the signature by feature matrix (basis) to
be greater than that which could be obtained by randomly shuffling the
input matrix.

    Step 5. Generate the chosen number of signatures with NMF. This is a heavy step and should be run in R instead of markdown mode. The output NMF contains sample by component matrix and component by signature matrix.

    Step 6: Quantify signatures: Given a sample-by-component matrix this function quantifies signature exposures using the LCD function from the YAPSA package, returning a normalized signature-by-sample matrix. If the component_by_signature matrix is specified then this matrix is used to define the signatures otherwise the signature definitions from the manuscript (feat_sig_mat.rds) are used. Signature exposure matrices were normalised to sum to one and exposures less than 0.01 were considered 0.

#### 7. CNsignature comparison between BRITROC, our HGSC tumor and HGSC Cervical samples:

I validate the TuCNsig which generated using Brenton's 36 CN_components
on either HGSC tumor samples (filtered away samples have zero tumor
fraction, 14 tumor samples) or HGSC Cervical samples (VS):

	1.  Generate CNsig from the HGSC tumor samples.
	2.  Generate CNsig from HGSC Cervical samples (VS).
	3.  Reorder component and plot heatmap.
	4.  Compare CNsignature across Britroc, HGSC tumor and HGSC Cervical.
	5.  Compare CNsig component weights: one histogram for each signature, containing the relative weighting of each of the components, colour coded by the feature distribution they come from.
	6.  Underlying Feature distributions of different sample sets.
	7.  Mixture model
	8.  CN Component means and standard deviations.
	10. Normalized signature exposure across cohorts
Here we compare the average exposure of each signature from different data sets or cohorts. 

#### 8. [CerCNsig_all](https://github.com/NyKepler/CerCNsig/tree/main/Tools/CNsignatures/Archived) validation and comparison Version 1

Next we will perform similar validation for the CerCNsig. The CerCNsig version 1 using all HGSC Cervical samples to extract the components (32) and signatures (6).
	
 	1.  Generate CerCNsig_all from the HGSC tumor samples.
	2.  Generate CerCNsig from HGSC tumor samples contains TF.
	3.  Generate CerCNsig from BritROC hq samples.
	4.  Reorder component and plot heatmap.
	5.  Compare CerCNsig_all across Britroc, HGSC tumor and HGSC Cervical
	6.  Compare CerCNsig_all component weights: one histogram for each signature, containing the relative weighting of each of the components, colour coded by the feature distribution they come from.
	7.  Compare CerCNsig_all to the CNsig
	8.  Underlying Feature distributions of different sample sets
	9.  Mixture model for CerCNsig_all
	10. VS-CN Component means and standard deviations
	11. Normalized signature exposure across cohorts

#### 9. [CerCNsig_filt](https://github.com/NyKepler/CerCNsig/tree/main/Tools/CNsignatures) validation and comparison Latest Version 2024 (ongoing)
The CerCNsig_filt version 2024 based on the absolution copy number profiles of HGSC Cervical samples generated from the https://github.com/IngridHLab/BINP52_CNA_Framework pipeline. Cervical samples were selected based on their HGSC CN signatures in Macintyre et al. 2018 https://github.com/markowetzlab/CNsignatures: samples with similiarity more than the first three signatures (S1-S3). Those cervical samples were considered to be CNA enriched instead of filtering the cervical samples using the cellularity from ACE/Rascal/ichorCNA estimation and mauanlly inspection which could be not completely accurate.
	
 	1. Filter away noisy cervical samples based on QDNAseq CN profile and the differences between expected standard deviation (Eσ) and measured standard deviation. 
  	2. Select the cervical samples have HGSC CNsig > 3 (42 samples remained).
   	3. 
   
  

#### 10. Assign CerCNsig on Benigh or BRCA Cervical samples.

#### 11. [Tools/BICseq2](https://github.com/NyKepler/CerCNsig/tree/main/Tools/BICseq2) copy number analysis: 
-    Most accurate annotation and copy number calling.

*BICseq2* is an algorithm developed for the normalization of
high-throughput sequencing (HTS) data and detect copy number variations
(CNV) in the genome. BICseq2 can be used for detecting CNVs with or
without a control genome. There are two main components in the
algorithm:

##### 1. Normalizing potential biases in the sequencing data using BICseq2-norm function. 

According to the BICseq2 requirements, I adjusted the *BWA* alignment
setting and *samtools* view function in a slightly different way in
order to filter out unique mapped read pairs with read length not
smaller than 100bp. *samtools* stats function provides the estimated
fragment size for later normalization.

    After aligment: Filter unique mapped read-pairs with same read length at least 100bp
    samtools view -@ $threads -F 260 -f 3 -q 1 -h $swd/$library.sort.bam | awk '/^@/ || length($10) >= 100' | samtools view -@ $threads -Sb > $swd/$library.filt.bam

    Prepocess: Generate seq file
    samtools view -@ $threads $bam/$library.filt.bam chr$x | perl $samtool_gu unique - | cut -f 4 > $readPosFile
    
The .seq files of each sample were then process by BICseq2-norm.pl to
remove the GC and mappability biases in the reads, and converted into
.txt files which contain 'obs' observed and 'expected' number of reads,
'var' variance and 'gc' GC percentage in every initial 10 bp bins along
each chromosome. Unlike other algorithms, BIC-seq2 performs
normalization at a nucleotide level (every 10 base-pairs) rather than at
a large bin level, resulting in its high sensitivity of detection for
small CNVs. Those inital bins will be used to calculate the Bayesian
information criterion in the segmentation steg.

```{bash BICseq2-norm}
## Options:
  	# --help
  	# -l=<int>: read length
  	# -s=<int>: fragment size
  	# -p=<float>: a subsample percentage: default 0.0002.
  	# -b=<int>: bin the expected and observed as <int> bp bins; Default 100.
  	# --gc_bin: if specified, report the GC-content in the bins
  	# --NoMapBin: if specified, do NOT bin the reads according to the mappability
  	# --bin_only: only bin the reads without normalization
	# --fig=<string>: plot the read count VS GC figure in the specified file (in pdf format)
  	# --title=<string>: title of the figure
  	# --tmp=<string>: the tmp directory;

## Generate normalized bin files
perl $BICSEQ_NORM -b=$bin -l=$readlen -s=$fragsize --gc_bin --NoMapBin --fig=$library --title=$library --tmp=$tmp $config $out
```
##### 2. Detecting CNVs based on the normalized data by using BICseq2-seg.

The segmentation can either run on stand-alone tumor sample or together
with a control sample ex. normal tissue or blood from the same patient.
Below I demonstrate the process with the control sample.

```{bash BICseq2-Segment}
## Options:
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

## For stand-alone sample without control!!
perl $BICSEQ_SEG --tmp=$tmp --lambda=$lambda --detail --bootstrap --fig=$png --title=$library.$binsize.$lambda $config $out 
```




