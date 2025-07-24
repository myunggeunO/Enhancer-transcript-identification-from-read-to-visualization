#!/bin/bash
set -e
THREADS=5

# Base directory setup
base_dir="01.aggregation_and_heatmap"
inter_dir="$base_dir/02.bed/inter"
intra_dir="$base_dir/02.bed/intra"
bw_dir="$base_dir/01.bigwig"

# 1. Heatmap
echo "Create Heatmap"
for subset in ES_active_E_intergenic ES_E_intergenic ES_non_active_E_intergenic Strand_ES_active_E_intragenic Strand_ES_E_intragenic Strand_ES_non_active_E_intragenic ; do
    if [[ "$subset" == *intergenic* ]]; then
        summit_bed="$inter_dir/${subset}_summit.bed"
        zMax_GRO=300
        zMax_H3K27ac=20
    else
        summit_bed="$intra_dir/${subset}_summit.bed"
        zMax_GRO=80
        zMax_H3K27ac=15
    fi

    matrix_out="$base_dir/03.Matrix/${subset}_heatmap_matrix.gz"
    GRO_out="$base_dir/04_2.heatmap/${subset}_GRO_heatmap.pdf"
    sorted_bed="$base_dir/02.bed/${subset}_descend.bed"
    sorted_out="$base_dir/03.Matrix/${subset}_sorted_matrix.gz"
    H3K27ac_out="$base_dir/04_2.heatmap/${subset}_H3K27ac_heatmap.pdf"

    # Matrix of GRO-seq generation
    computeMatrix reference-point --referencePoint center -a 1000 -b 1000 -R "$summit_bed" -S "$bw_dir/GRO_inter_intraE_signal.bw" --missingDataAsZero -p "$THREADS" -o "$matrix_out"

    # Heatmap of GRO-seq generation
    plotHeatmap -m "$matrix_out" -out "$GRO_out" --dpi 500 --sortRegions descend --refPointLabel "summit" --legendLocation best --samplesLabel "GRO-seq" --whatToShow "heatmap and colorbar" --heatmapHeight 10 --heatmapWidth 6 --colorList "w,#f5a118" --zMin 0 --zMax "$zMax_GRO" --outFileSortedRegions "$sorted_bed"

    # Matrix of H3K27ac generation sorted by GRO signal
    computeMatrix reference-point --referencePoint center -a 1000 -b 1000 -R "$sorted_bed" -S "$bw_dir/H3K27ac.bw" --missingDataAsZero -p "$THREADS" -o "$sorted_out"

    # Heatmap of H3K27ac generation sorted by GRO signal
    plotHeatmap -m "$sorted_out" -out "$H3K27ac_out" --dpi 500 --sortRegions keep --refPointLabel "summit" --legendLocation best --samplesLabel "H3K27ac" --whatToShow "heatmap and colorbar" --colorList "w,#e71e15" --heatmapHeight 10 --heatmapWidth 6 --zMin 0.8 --zMax "$zMax_H3K27ac"

done

echo ""
echo "Generation of heatmap done"
