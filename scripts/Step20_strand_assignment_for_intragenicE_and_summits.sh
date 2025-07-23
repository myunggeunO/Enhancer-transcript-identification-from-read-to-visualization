#!/bin/bash
set -e

# Define base paths
input_all="03.intraE/strand_designation/01.overlapped_gene_strand/ES_E_intragenic_strand_with_dup.bed"
strand_bed="03.intraE/strand_designation/04.enhancer_strand_designation_by_RPKM_of_gene/ES_E_intragenic_PCG_priority_strand_by_gene_RPKM.bed"

# Define output directories
out_dir="03.intraE/final_strand_IntragenicE"
summit_dir="03.intraE/summit"
strand_summit_dir="$out_dir/summit"
mkdir -p "$strand_summit_dir"

# 1. Subtract opposite strand from all strand-assigned enhancers
echo "Subtract opposite-strand regions from strand unfiltered enhancers"
bedtools subtract -a "$input_all" -b "$strand_bed" -S > "$out_dir/Strand_ES_E_intragenic.bed"

# 2. Assign strand to active/non-active enhancer sets
echo "Assign strand to active intragenic enhancers"
bedtools intersect -wa -u -a "$out_dir/Strand_ES_E_intragenic.bed" -b "03.intraE/ES_active_E_intragenic.bed" > "$out_dir/Strand_ES_active_E_intragenic.bed"

echo "Assign strand to non-active intragenic enhancers"
bedtools intersect -wa -u -a "$out_dir/Strand_ES_E_intragenic.bed" -b "03.intraE/ES_non_active_E_intragenic.bed" > "$out_dir/Strand_ES_non_active_E_intragenic.bed"

# 3. Assign strand to enhancer summit files
echo "Assign strand to summit files for each enhancer set"
for type in ES_E_intragenic ES_active_E_intragenic ES_non_active_E_intragenic; do
    summit_in="$summit_dir/${type}_summit.bed"
    summit_out="$strand_summit_dir/Strand_${type}_summit.bed"

    bedtools intersect -wa -wb -a "$summit_in" -b "$out_dir/Strand_ES_E_intragenic.bed" | \
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$NF}' > "$summit_out"
done

echo ""
echo "Strand assignment for all enhancer sets and summit files done"
