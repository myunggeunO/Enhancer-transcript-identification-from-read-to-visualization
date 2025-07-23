#!/bin/bash
set -e

# Set base directory
base_dir="01.E_identification/material"
mkdir -p "$base_dir"/{ATAC,Histone,Annotation}

# Copy ATAC peaks
echo "Copy ATAC peaks"
cp MATERIAL/ATAC/peak_calling/ATAC_peaks.narrowPeak "$base_dir/ATAC/"
cp MATERIAL/ATAC/peak_calling/ATAC_summits.bed "$base_dir/ATAC/"

# Copy histone peaks
echo "Copy histone peaks"
cp MATERIAL/H3K4me1/peak_calling/overlapped_peak/H3K4me1_overlapped_peaks.broadpeak "$base_dir/Histone/"
cp MATERIAL/H3K27ac/peak_calling/H3K27ac_peaks.broadPeak "$base_dir/Histone/"

# Download GENCODE annotation
echo "Download GENCODE annotation"
wget -c -P "$base_dir/Annotation" \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz

# Copy Blacklist of mm10
cp mm10-blacklist/mm10-blacklist.v2.bed "$base_dir/ATAC/"
cp mm10.chromsizes/mm10.chrom.sizes "$base_dir/Annotation"

echo ""
echo "Ready for enhancer identification"
