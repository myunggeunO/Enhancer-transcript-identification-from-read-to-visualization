#!/bin/bash
set -e
THREADS=5

# Set directories
region_dir="02.eRNA_quantification/02.region"
bam="02.eRNA_quantification/01.bam/SRR5655667_filtered.sorted.bam"
out_root="02.eRNA_quantification/03.count_normalized_with_RPKM"
mkdir -p "$out_root"/{inter,intra}

# 1. Output total mapped reads of GRO-seq
echo "Extract total mapped reads of GRO-seq"
total=$(sambamba flagstat "$bam" | awk '/in total/ {print $1}')
echo "$total" > "$out_root/GRO_total_reads.txt"

# 2. Convert BED to GTF
echo "Convert BED file to GTF format for featureCounts"
for bed in "$region_dir"/{inter,intra}/*.bed; do
    base=$(basename "$bed" .bed)
    gtf="${bed%.bed}.gtf"
    awk -F "\t" 'BEGIN{OFS="\t"}{print $1,"bed","enhancer",$2+1,$3,".",$6,".","gene_id \""$4"\";"}' "$bed" > "$gtf"

# 3. Count mapped reads of enhancer sets
if [[ "$bed" == *"/intra/"* ]]; then
    strand=2; out_dir="$out_root/intra"; name="${base}_AS"
else
    strand=0; out_dir="$out_root/inter"; name="${base}"
fi
echo "Run featureCounts to count mapped reads of $name"
counts="${out_dir}/${name}_counts.txt"
rpkm_out="${out_dir}/${name}_RPKM.txt"
featureCounts -a "$gtf" -o "$counts" -T "$THREADS" -s "$strand" -t enhancer -g gene_id -O --fraction "$bam"

# 4. Calculte RPKM of enhancer mapped reads
awk -v total="$total" 'BEGIN{OFS="\t"} NR>2 {
  len=$6; count=$7; kb=len/1000
  if (kb>0 && total>0) {
    rpk=count/kb; rpkm=rpk/(total/1e6)
    print $1, count, len, rpk, rpkm
  } else {
    print $1, count, "N/A", "N/A", "N/A"
  }
}' "$counts" > "$rpkm_out"
done

echo ""
echo "All RPKM calculations completed and organized by inter/intra"
