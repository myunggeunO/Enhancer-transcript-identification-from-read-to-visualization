#!/bin/bash
set -e
THREADS=5

# Define base paths
base_dir="03.intraE/strand_designation"
rpkm_dir="$base_dir/03.RPKM_calculation_of_overlapped_gene"
mkdir -p "$rpkm_dir"

# Input files
enhancer_bed="$base_dir/02.enhancer_with_PCG_priority/ES_E_intragenic_PCG_priority.bed"
gene_bed="material/Annotation/M23_gene.bed"
bam="../MATERIAL/GRO/rep1/02.Align/SRR5655667_filtered.sorted.bam"

# Output files
overlap_bed="$rpkm_dir/overlapped_genes.bed"
gtf="$rpkm_dir/overlapped_genes.gtf"
counts="$rpkm_dir/overlapped_gene_counts.tsv"
total_txt="$rpkm_dir/GRO_total_reads.txt"
rpkm_out="$rpkm_dir/overlapped_gene_RPKM.tsv"

# 1. Find strand-specific overlapping genes
echo "Print gene list overlapping with initially strand-assigned enhancers"
bedtools intersect -s -wa -u -a "$gene_bed" -b "$enhancer_bed" > "$overlap_bed"

# 2. Convert BED to GTF
echo "Convert BED file to GTF format for featureCounts"
awk 'BEGIN{OFS="\t"}{print $1,"bed","gene",$2+1,$3,".",$6,".","gene_id \""$4"\";"}' "$overlap_bed" > "$gtf"

# 3. Count reads
echo "Run featureCounts to count mapped reads per gene"
featureCounts -a "$gtf" -o "$counts" -T "$THREADS" -s 1 -t gene -g gene_id -O --fraction "$bam"

# 4. Get total reads
sambamba flagstat "$bam" | awk '/in total/ {print $1}' > "$total_txt"
total=$(cat "$total_txt")

# 5. Calculate RPKM
echo "Compute RPKM from GTF, mapped counts, and total read counts"
awk -v total="$total" 'BEGIN{OFS="\t"} NR>2 {
  len=$6; count=$7; kb=len/1000
  if (kb>0 && total>0) {
    rpk=count/kb; rpkm=rpk/(total/1e6)
    print $1, count, len, rpk, rpkm
  } else {
    print $1, count, "N/A", "N/A", "N/A"
  }
}' "$counts" > "$rpkm_out"

echo ""
echo "RPKM calculation for second strand assignment done"
