#!/bin/bash
set -e

# 1. Detect Mamba/Conda
if command -v mamba >/dev/null 2>&1; then
  CREATE_CMD="mamba"
elif command -v conda >/dev/null 2>&1; then
  CREATE_CMD="conda"
else
  echo "ERROR: Neither 'mamba' nor 'conda' found in PATH. Follow Step1."
  exit 1
fi

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

# 2. Install bioinformatics tools
${CREATE_CMD} install -c bioconda -c conda-forge sra-tools && echo "sra-tools (fasterq-dump, prefetch) installed"
${CREATE_CMD} install -c bioconda -c conda-forge trim-galore && echo "trim-galore installed"
${CREATE_CMD} install -c conda-forge pigz && echo "pigz installed"
${CREATE_CMD} install -c bioconda -c conda-forge bowtie2 -y && echo "bowtie2 installed"
${CREATE_CMD} install -c bioconda -c conda-forge bedtools -y && echo "bedtools installed"
${CREATE_CMD} install -c bioconda -c conda-forge cutadapt -y && echo "cutadapt installed"
${CREATE_CMD} install -c bioconda -c conda-forge fastqc -y && echo "fastqc installed"
${CREATE_CMD} install -c bioconda -c conda-forge macs3 -y && echo "macs3 installed"
${CREATE_CMD} install -c bioconda -c conda-forge subread -y && echo "featureCounts (subread) installed"
${CREATE_CMD} install -c bioconda -c conda-forge samtools -y && echo "samtools installed"
${CREATE_CMD} install -c bioconda -c conda-forge sambamba -y && echo "sambamba installed"
${CREATE_CMD} install -c bioconda -c conda-forge deeptools -y && echo "deeptools installed"
${CREATE_CMD} install -c bioconda -c conda-forge ucsc-bedgraphtobigwig -y && echo "bigWigbuilder (bedgraphtobigwig) installed"

# 2. Install R packages
${CREATE_CMD} install -c conda-forge r-tidyverse -y && echo "r-tidyverse installed"
${CREATE_CMD} install -c conda-forge r-cowplot -y && echo "r-cowplot installed"

# 3. Install HOMER
echo "!Installing HOMER..."
mkdir -p ~/homer && cd ~/homer
wget -q http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install
perl configureHomer.pl -install mm10 && echo "HOMER installed"

echo ""
echo "All packages installed successfully in enhancer-env"
