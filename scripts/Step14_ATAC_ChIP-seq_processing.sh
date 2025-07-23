#!/bin/bash
set -e

# Define directories
base_dir="material"
atac_dir="$base_dir/ATAC"
histone_dir="$base_dir/Histone"
anno_dir="$base_dir/Annotation"
promoter="$anno_dir/M23_promoter_candidate.bed"
chromsizes="$anno_dir/mm10.chrom.sizes"
blacklist="$atac_dir/mm10-blacklist.v2.bed"

# 1. Remove blacklist and promoters from ATAC

echo "Process ATAC peaks by removing blacklist and promoter"

for bed in "$atac_dir"/*narrowPeak; do
    [[ ! -f $bed ]] && continue
    mark="ATAC"

    out_blacklist="$atac_dir/${mark}_no_blacklist.bed"
    bedtools subtract -A -a "$bed" -b "$blacklist" > "$out_blacklist"

    echo "Remove promoters from $out_blacklist"
    out_promoter="$atac_dir/${mark}_no_blacklist_no_promoter.bed"
    bedtools subtract -A -a "$out_blacklist" -b "$promoter" > "$out_promoter"
done

# 2. Expand histone peaks ±1 kb and exclude promoter regions

echo "Process histone peaks with +-1kb slop and promoter filtering"

flanking_dir="$histone_dir/1kb_flanking"
filtered_dir="$histone_dir/1kb_flanking_no_promoter"
mkdir -p "$flanking_dir" "$filtered_dir"

for mark in H3K4me1 H3K27ac; do
    echo "Process $mark"

    input_file="$histone_dir/${mark}_overlapped_peaks.broadpeak"
    [[ "$mark" == "H3K27ac" ]] && input_file="$histone_dir/${mark}_peaks.broadPeak"

    out_slop="$flanking_dir/${mark}_+-1kb.bed"
    out_final="$filtered_dir/${mark}_+-1kb_no_promoter.bed"

    echo "Apply ±1kb flanking"
    bedtools slop -b 1000 -g "$chromsizes" -i "$input_file" > "$out_slop"

    echo "Remove promoters"
    bedtools subtract -a "$out_slop" -b "$promoter" > "$out_final"
done

echo ""
echo "All processing complete."
