#!/bin/bash
set -e

# Define input files
intragenic_enhancer="03.intraE/ES_E_intragenic.bed"
gene_annotation="material/Annotation/M23_gene.bed"

# Define output directories
out_dir="03.intraE/strand_designation/01.overlapped_gene_strand"
mkdir -p "$out_dir"

# 1. Intersect enhancer with gene body
echo "Find overlaps between intragenic enhancers and gene annotations"
bedtools intersect -wa -wb -a "$intragenic_enhancer" -b "$gene_annotation" > "${out_dir}/enhancer_gene_pair.tsv"

# 2. Extract specific fields and remove duplicates
# Fields: enhancer_chr, start, end, id, enhancer_strand, gene_strand
echo "Extract relevant columns and removing duplicates"
awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $16}' "${out_dir}/enhancer_gene_pair.tsv" | sort -k1,1 -k2,2n | uniq > "${out_dir}/ES_E_intragenic_strand_with_dup.bed"

echo ""
echo "Tempotal strand designation complete."
