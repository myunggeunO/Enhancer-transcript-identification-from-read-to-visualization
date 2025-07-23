#!/bin/bash
set -e

# 1. Create reference index directory
mkdir -p reference_index
cd reference_index

# 2. Download Bowtie2 index for mm10
wget https://genome-idx.s3.amazonaws.com/bt/mm10.zip

# 3. Unzip the index archive
unzip mm10.zip

# 4. Remove the zip file
rm mm10.zip

# 5. Print completion message
echo "mm10 Bowtie2 index ready"
