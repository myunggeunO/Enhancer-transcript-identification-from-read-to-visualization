#!/bin/bash
set -e

# Define input files
atac_peaks="material/ATAC/ATAC_no_blacklist_no_promoter.bed"
atac_summits="material/ATAC/ATAC_summits.bed"
k4me1="material/Histone/1kb_flanking_no_promoter/H3K4me1_+-1kb_no_promoter.bed"
k27ac="material/Histone/1kb_flanking_no_promoter/H3K27ac_+-1kb_no_promoter.bed"
gene_body="material/Annotation/M23_gene.bed"

# Define output directories
mkdir -p 01.allE/summit 02.interE/summit 03.intraE/summit

# 1. Identify ES_E (ATAC ∩ H3K4me1)
echo "Identify ES_E (ATAC ∩ H3K4me1)"
bedtools intersect -wa -u -a "$atac_peaks" -b "$k4me1" > 01.allE/ES_E.bed

# 2. Identify ES_active_E (ES_E ∩ H3K27ac)
echo "Identify ES_active_E (ES_E ∩ H3K27ac)"
bedtools intersect -wa -u -a 01.allE/ES_E.bed -b "$k27ac" > 01.allE/ES_active_E.bed

# 3. Identify ES_non_active_E (ES_E - ES_active_E)
echo "Identify ES_non_active_E (ES_E - ES_active_E)"
bedtools subtract -a 01.allE/ES_E.bed -b 01.allE/ES_active_E.bed > 01.allE/ES_non_active_E.bed

# 4. Count validation
total=$(wc -l < 01.allE/ES_E.bed)
active=$(wc -l < 01.allE/ES_active_E.bed)
nonactive=$(wc -l < 01.allE/ES_non_active_E.bed)

echo "!!Enhancer counts → Total: $total = Active: $active + Non-active: $nonactive"
if [[ $((active + nonactive)) -ne $total ]]; then
    echo "[ERROR] Enhancer count mismatch"
    exit 1
fi

# 5. Extract summit for each enhancer type from ATAC_summits
for enh in ES_E ES_active_E ES_non_active_E; do
    echo "Extract summit for $enh"
    bedtools intersect -u -a "$atac_summits" -b "01.allE/${enh}.bed" > "01.allE/summit/${enh}_summit.bed"
done

# 6. Split each summit into intergenic/intragenic using gene body
for summit in 01.allE/summit/*_summit.bed; do
    base=$(basename "$summit" .bed)     # e.g., ES_E_summit
    name=${base%_summit}                # e.g., ES_E

    echo "Split summit: $name"
    bedtools intersect -v -a "$summit" -b "$gene_body" > "02.interE/summit/${name}_intergenic_summit.bed"
    bedtools intersect -u -a "$summit" -b "$gene_body" > "03.intraE/summit/${name}_intragenic_summit.bed"
done

# 7. Map enhancer regions to inter/intra based on summit location
for name in ES_E ES_active_E ES_non_active_E; do
    echo "Map enhancer by summit: $name"

    inter_summit="02.interE/summit/${name}_intergenic_summit.bed"
    intra_summit="03.intraE/summit/${name}_intragenic_summit.bed"

    bedtools intersect -u -a "01.allE/${name}.bed" -b "$inter_summit" > "02.interE/${name}_intergenic.bed"
    bedtools intersect -u -a "01.allE/${name}.bed" -b "$intra_summit" > "03.intraE/${name}_intragenic.bed"
done

echo ""
echo "Enhancer identification and classification complete."
