#!/bin/bash
set -e
THREADS=5

# Base directory setup
base_dir="01.aggregation_and_heatmap"
inter_dir="$base_dir/02.bed/inter"
intra_dir="$base_dir/02.bed/intra"
bw_dir="$base_dir/01.bigwig"

# 1. aggregation plot
echo "Create aggregation plot"
for subset in ES_active_E_intergenic ES_E_intergenic ES_non_active_E_intergenic Strand_ES_active_E_intragenic Strand_ES_E_intragenic Strand_ES_non_active_E_intragenic; do
    if [[ "$subset" == *intergenic* ]]; then
        summit_bed="$inter_dir/${subset}_summit.bed"
    else
        summit_bed="$intra_dir/${subset}_summit.bed"
    fi
    matrix_out="$base_dir/03.Matrix/${subset}_aggregation_matrix.gz"
    profile_out="$base_dir/04_1.aggregation/${subset}_aggregation.pdf"

    # Plot generagion
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$summit_bed" \
        -S "$bw_dir/ATAC.bw" "$bw_dir/H3K27ac.bw" "$bw_dir/H3K4me1.bw" \
        --missingDataAsZero \
        -p "$THREADS" \
        -o "$matrix_out"

    plotProfile \
        -m "$matrix_out" \
        -out "$profile_out" \
        --plotHeight 8 \
        --plotWidth 8 \
        --yMax 6 \
        --perGroup \
        --refPointLabel "summit" \
        --regionsLabel "${subset}" \
        --samplesLabel ATAC H3K27ac H3K4me1 \
        --colors '#2ca339' '#e7211a' '#f4da26'

done

echo ""
echo "Generation of aggregation plot complete"
