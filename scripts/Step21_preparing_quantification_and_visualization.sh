#!/bin/bash
set -e
THREADS=5

# Base directory setup
base_dir="02.E_visualization_quantification"
agg_dir="$base_dir/01.aggregation_and_heatmap"
bigwig_dir="$agg_dir/01.bigwig"
bed_dir="$agg_dir/02.bed"
matrix_dir="$agg_dir/03.Matrix"
agg_img_dir="$agg_dir/04_1.aggregation"
heatmap_dir="$agg_dir/04_2.heatmap"
erna_quant_dir="$base_dir/02.eRNA_quantification"
erna_bam_dir="$erna_quant_dir/01.bam"
erna_region_dir="$erna_quant_dir/02.region"
erna_count_dir="$erna_quant_dir/03.count_normalized_with_RPKM"
erna_vis_dir="$base_dir/03.eRNA_visualization"

# Make path for visualization and quantification
mkdir -p "$bigwig_dir" "$bed_dir/inter" "$bed_dir/intra" "$matrix_dir" "$agg_img_dir" "$heatmap_dir"
mkdir -p "$erna_bam_dir" "$erna_region_dir/inter" "$erna_region_dir/intra" "$erna_count_dir" "$erna_vis_dir"

# Copy summit files
inter_summit_src="01.E_identification/02.interE/summit"
intra_summit_src="01.E_identification/03.intraE/final_strand_IntragenicE/summit"

echo "Copy summit files"
cp "$inter_summit_src"/* "$bed_dir/inter/"
cp "$intra_summit_src"/* "$bed_dir/intra/"

# Copy bigwig files
atac_bw="MATERIAL/ATAC/merge/04.bigwig/ATAC.bw"
h3k27ac_bw="MATERIAL/H3K27ac/tmerge/04.bigwig/H3K27ac.bw"
h3k4me1_bw="MATERIAL/H3K4me1/merge/04.bigwig/H3K4me1.bw"

echo "Copy bigwig signals"
cp "$atac_bw" "$bigwig_dir/"
cp "$h3k27ac_bw" "$bigwig_dir/"
cp "$h3k4me1_bw" "$bigwig_dir/"

# Copy GRO-seq bam and enhancers to eRNA quantification path
gro_bam="MATERIAL/GRO/rep1/02.Align/SRR5655667_filtered.sorted.bam"

echo "Copy enhancer files and GRO-seq bam file to eRNA quantification path"
cp "$inter_summit_src/../"*_active_*.bed "$erna_region_dir/inter/"
cp "$intra_summit_src/../"*_active_*.bed "$erna_region_dir/intra/"
cp $gro_bam $erna_bam_dir/

echo ""
echo "Quantification and visualization materials are ready"
