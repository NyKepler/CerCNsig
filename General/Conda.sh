#!/usr/bin/env bash

# Install Miniconda3 
# mkdir /home/minerva/Downloads
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /home/minerva/Downloads/Miniconda.sh
chmod +x /home/minerva/Downloads/Miniconda.sh
/home/minerva/Downloads/Miniconda.sh #choose directory ex. /Users/minerva/miniconda3
source /home/minerva/miniconda3/bin/activate
conda init

# Setup channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Create enviroments and install bwa, samtools, Picardtools and gatk
conda create -n sWGS1.0 bwa samtools gatk4 picard trim-galore 
conda activate sWGS1.0 
conda install sra-tools=2.9
# vdb-config -i
conda install python multiqc
conda install -c hcc gistic2 
conda install -c bioconda bioconductor-cntools bedtools
conda install -c bioconda hmmcopy #for ichorcna only generate read counts (same as HMMcopy suite)
conda install -c bioconda freebayes #VCF calling
conda install -c bioconda bcftools #VCF calling
# mutect2 is installed through gatk4 package
conda install -c bioconda strelka #VCF calling, 2.9.10 version
conda install -c bioconda nextflow #sarek pipeline.
conda deactivate


#Install R
conda create -n R r-essentials #r-rgraphics can't be installed
conda activate R
conda install r-r.cache r-digest 
conda install bioconductor-qdnaseq bioconductor-qdnaseq.hg19 bioconductor-dnacopy bioconductor-cghcall bioconductor-cghregions #
conda install r-nmf r-flexmix bioconductor-yapsa r-domc #CNSignature
conda install -c bioconda bioconductor-ace #ACE for absolute cna
conda install -c bioconda r-ichorcna bioconductor-hmmcopy #ichorcna (latest version) and hmmcopy for CN analysis.
conda install -c bioconda bioconductor-bsgenome.hsapiens.ucsc.hg19 ucsc-bigwigaverageoverbed  #QDNAseq bin annotation
conda install -c conda-forge r-devtools r-stringi r-stringr #Rascal for absolute cna
conda install -c bioconda bioconductor-karyoploter bioconductor-regioner #genomoes displaying arbitrary data
conda install -c conda-forge r-ggpubr

 
# nf-core 
conda create --name nf-core python=3.8 nf-core nextflow

