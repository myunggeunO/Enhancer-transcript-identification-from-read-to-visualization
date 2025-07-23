#!/bin/bash
set -e

# Define key paths
base_dir="03.intraE/strand_designation"
input_bed="$base_dir/01.overlapped_gene_strand/ES_E_intragenic_strand_with_dup.bed"
pcg_bed="material/Annotation/M23_protein_coding_gene.bed"

# Output setup
out_dir="$base_dir/02.enhancer_with_PCG_priority"
mkdir -p "$out_dir"
prefix="$out_dir/ES_E_intragenic"

# 1. Select enhancers present on both + and - strands
echo "Filter enhancers present on both strands"
awk '
{
    key = $1"\t"$2"\t"$3"\t"$4"\t"$5
    strand = $6
    data[key][strand] = (key in data && strand in data[key]) ? data[key][strand] ORS $0 : $0
    strand_seen[key][strand] = 1
}
END {
    for (k in strand_seen)
        if ("+" in strand_seen[k] && "-" in strand_seen[k]) {
            print data[k]["+"]
            print data[k]["-"]
        }
}
' "$input_bed" > "${prefix}_both_strands.bed"

# 2. Same-strand overlap with PCG
echo "Identify PCG-overlapping enhancers depending on strand"
bedtools intersect -s -wa -u -a "${prefix}_both_strands.bed" -b "$pcg_bed" > "${prefix}_overlap_PCG_same_strand.bed"

# 3. No overlap with PCG
echo "Identify enhancers not overlapping any PCG"
bedtools intersect -v -a "${prefix}_both_strands.bed" -b "$pcg_bed" > "${prefix}_no_PCG_overlap.bed"

# 4. Combine and sort
echo "Merge PCG-prioritized and non-overlapping enhancers"
cat "${prefix}_overlap_PCG_same_strand.bed" "${prefix}_no_PCG_overlap.bed" | sort -k1,1 -k2,2n > "${prefix}_PCG_priority.bed"

echo ""
echo "Strand assignment for enhancers overlapping genes on both strands with PCG prioritization done"
