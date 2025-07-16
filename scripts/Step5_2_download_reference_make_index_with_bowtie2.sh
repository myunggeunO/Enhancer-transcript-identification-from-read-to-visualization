#!/bin/bash
set -e

# Set the number of threads
THREADS=5

# 1. Create reference index directory
mkdir -p reference_index
cd reference_index

# 2. Download mm10 FASTA file from UCSC
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz

# 3. Decompress the FASTA file using pigz
pigz -p "$THREADS" -d mm10.fa.gz

# 4. Build Bowtie2 index from FASTA
bowtie2-build mm10.fa mm10

# 5. Print completion message
echo "Bowtie2 index built from UCSC mm10.fa"
