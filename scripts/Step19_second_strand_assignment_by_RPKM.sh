#!/bin/bash
set -e

# Define base paths
base_dir="03.intraE/strand_designation"
rpkm_dir="$base_dir/04.enhancer_strand_designation_by_RPKM_of_gene"
mkdir -p "$rpkm_dir"

# Input files
enhancer_bed="$base_dir/02.enhancer_with_PCG_priority/ES_E_intragenic_PCG_priority.bed"
gene_bed="$base_dir/03.RPKM_calculation_of_overlapped_gene/overlapped_genes.bed"
rpkm_file="$base_dir/03.RPKM_calculation_of_overlapped_gene/overlapped_gene_RPKM.tsv"

# Output files
gene_rpkm_bed="$rpkm_dir/gene_with_rpkm.bed"
merged="$rpkm_dir/enhancer_gene_merged.bed"
final_bed="$rpkm_dir/ES_E_intragenic_PCG_priority_strand_by_gene_RPKM.bed"

# 1. Join gene BED with RPKM info
echo "Join gene annotation with RPKM values"
awk '{split($1,a,"."); print a[1]"\t"$5}' "$rpkm_file" | sort -k1,1 > "$rpkm_dir/sorted_rpkm.txt"
awk 'BEGIN{OFS="\t"} {split($4,a,"."); print a[1], $0}' "$gene_bed" | sort -k1,1 > "$rpkm_dir/sorted_genes.txt"
join -1 1 -2 1 -t $'\t' "$rpkm_dir/sorted_rpkm.txt" "$rpkm_dir/sorted_genes.txt" | \
awk 'BEGIN{OFS="\t"}{print $3,$4,$5,$6,$2,$8}' > "$gene_rpkm_bed"

# 2. Intersect enhancer with gene with RPKM info
echo "Intersect enhancers with gene annotation and RPKM values"
bedtools intersect -s -wa -wb -a "$enhancer_bed" -b "$gene_rpkm_bed" > "$merged"

# 3. Pick gene with highest RPKM per enhancer
echo "Select the overlapping gene with highest RPKM per enhancer"
awk -F'\t' '
{
    key = $1"\t"$2"\t"$3
    rpkm = $(NF-1) + 0
    if (!(key in best) || rpkm > best[key]) {
        best[key] = rpkm
        line[key] = $0
    }
}
END {
    for (k in line) print line[k]
}' "$merged" | \
awk 'BEGIN{OFS="\t"} {rpkm=$(NF-2); strand=$(NF); print $1,$2,$3,$4,rpkm,strand}' | sort -k1,1 -k2,2n > "$final_bed"

echo ""
echo "Enhancer strand assignment by gene RPKM done"
