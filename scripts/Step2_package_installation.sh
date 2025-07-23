#!/bin/bash

set -e

# List of bioinformatics tools to be installed
# - sra_tools
# - trim-galore
# - pigz
# - bowtie2
# - bedtools
# - cutadapt
# - fastqc
# - macs3
# - featureCounts (via subread)
# - samtools
# - sambamba
# - deeptools

# R packages:
# - r-tidyverse
# - r-cowplot

# Other:
# - HOMER

# 1. Install bioinformatics tools
mamba install -c bioconda -c conda-forge sra-tools && echo "sra-tools (fasterq-dump, prefetch) installed"
mamba install -c bioconda -c conda-forge trim-galore && echo "trim-galore installed"
mamba install -c conda-forge pigz && echo "pigz installed"
mamba install -c bioconda -c conda-forge bowtie2 -y && echo "bowtie2 installed"
mamba install -c bioconda -c conda-forge bedtools -y && echo "bedtools installed"
mamba install -c bioconda -c conda-forge cutadapt -y && echo "cutadapt installed"
mamba install -c bioconda -c conda-forge fastqc -y && echo "fastqc installed"
mamba install -c bioconda -c conda-forge macs3 -y && echo "macs3 installed"
mamba install -c bioconda -c conda-forge subread -y && echo "featureCounts (subread) installed"
mamba install -c bioconda -c conda-forge samtools -y && echo "samtools installed"
mamba install -c bioconda -c conda-forge sambamba -y && echo "sambamba installed"
mamba install -c bioconda -c conda-forge deeptools -y && echo "deeptools installed"

# 2. Install R packages
mamba install -c conda-forge r-tidyverse -y && echo "r-tidyverse installed"
mamba install -c conda-forge r-cowplot -y && echo "r-cowplot installed"

# 3. Install HOMER
echo "!Installing HOMER..."
mkdir -p ~/homer && cd ~/homer
wget -q http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install
perl configureHomer.pl -install mm10 && echo "HOMER installed"

echo ""
echo "All packages installed successfully in enhancer-env"
