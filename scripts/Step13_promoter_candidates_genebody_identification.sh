#!/bin/bash
set -e

#"material/annotation"
base_dir="material/Annotation"
gtf="${base_dir}/gencode.vM23.annotation.gtf"
chromsizes="${base_dir}/mm10.chrom.sizes"

# Decompress GTF if not already present
if [[ ! -f "$gtf" ]]; then
    echo "Decompressing GTF file"
    gunzip -c "$gtf.gz" > "$gtf"
fi

# Generate promoter candidate BED: TSS ±2kb (safe with bedtools slop)
echo "Generate M23_promoter_candidate.bed"
awk -F'\t' -v OFS='\t' '$3 == "transcript" {
    match($0, /gene_id "[^"]+"/, gid); gsub(/gene_id "|"/, "", gid[0]); gene_id = gid[0];
    if ($7 == "+") {tss = $4 - 1} else if ($7 == "-") {tss = $5} else {next}
    print $1, tss, tss + 1, gene_id, ".", $7
}' "$gtf" | bedtools slop -b 2000 -g "$chromsizes" | bedtools sort -g "$chromsizes" > "$base_dir/M23_promoter_candidate.bed"

echo "Promoter candidate done"

# Generate gene BED: full gene body from gene features
echo "Generate M23_gene.bed"
awk -F'\t' -v OFS='\t' '$3 == "gene" {
    match($0, /gene_id "[^"]+"/, gid); gsub(/gene_id "|"/, "", gid[0]); gene_id = gid[0];
    print $1, $4 - 1, $5, gene_id, ".", $7
}' "$gtf" | bedtools sort -g "$chromsizes" > "$base_dir/M23_gene.bed"

echo "Gene body done"

# Generate protein_coding gene BED
echo "Generate M23_protein_coding_gene.bed"
awk -F'\t' -v OFS='\t' '$3 == "gene" && $0 ~ /gene_type "protein_coding"/ {
    match($0, /gene_id "[^"]+"/, gid); gsub(/gene_id "|"/, "", gid[0]); gene_id = gid[0];
    print $1, $4 - 1, $5, gene_id, ".", $7
}' "$gtf" | bedtools sort -g "$chromsizes" > "$base_dir/M23_protein_coding_gene.bed"

echo "Protein coding gene body done"

echo ""
echo "All done"
