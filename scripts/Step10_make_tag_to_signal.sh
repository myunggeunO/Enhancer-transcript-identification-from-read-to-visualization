#!/bin/bash
set -e
CHROMSIZES="mm10.chromsizes/mm10.chrom.sizes"

# 1. GRO
echo "GRO_Tag_and_Signal"
bam=(MATERIAL/GRO/rep1/02.Align/*.sorted.bam)
tagdir="MATERIAL/GRO/rep1/03.TagDir/GRO"
[[ -f "${bam[0]}" ]] && {
    mkdir -p "$tagdir"
    makeTagDirectory "$tagdir" "${bam[0]}"

    outdir="MATERIAL/GRO/rep1/04.bigwig"
    mkdir -p "$outdir"
    for strand in + -; do
        makeUCSCfile "$tagdir" -style rnaseq -strand "$strand" -bigWig "$CHROMSIZES" -o "$outdir/GRO_${strand}.bw"
    done
} || echo "[WARNING] BAM not found for GRO (${bam[0]})"

# 2. H3K27ac input
echo "H3K27ac_Tag_and_Signal"
echo "Input_Tag"
bam="MATERIAL/H3K27ac/tmerge/02.Align/H3K27ac_input_merge_filtered_rmdup.sorted.bam"
input_tagdir="MATERIAL/H3K27ac/tmerge/03.TagDir/H3K27ac_input"
[[ -f "$bam" ]] && {
    mkdir -p "$input_tagdir"
    makeTagDirectory "$input_tagdir" "$bam"
} || echo "[WARNING] BAM not found for H3K27ac input ($bam)"

# 3. H3K27ac
echo "H3K27ac_Tag_Singnal"
bam="MATERIAL/H3K27ac/tmerge/02.Align/H3K27ac_merge_filtered_rmdup.sorted.bam"
tagdir="MATERIAL/H3K27ac/tmerge/03.TagDir/H3K27ac"
[[ -f "$bam" ]] && {
    mkdir -p "$tagdir"
    makeTagDirectory "$tagdir" "$bam"

    outdir="MATERIAL/H3K27ac/tmerge/04.bigwig"
    mkdir -p "$outdir"
    makeUCSCfile "$tagdir" -bigWig "$CHROMSIZES" -c "$input_tagdir" -pseudo 1 -o "$outdir/H3K27ac.bw"
} || echo "[WARNING] BAM not found for H3K27ac ($bam)"

# 4. ATAC
echo "ATAC_Tag_and_to_Signal"
bam="MATERIAL/ATAC/merge/02.Align/ATAC_merge.bam"
tagdir="MATERIAL/ATAC/merge/03.TagDir/ATAC"
[[ -f "$bam" ]] && {
    mkdir -p "$tagdir"
    makeTagDirectory "$tagdir" "$bam"

    outdir="MATERIAL/ATAC/merge/04.bigwig"
    mkdir -p "$outdir"
    makeUCSCfile "$tagdir" -bigWig "$CHROMSIZES" -o "$outdir/ATAC.bw"
} || echo "[WARNING] BAM not found for ATAC ($bam)"

# 5. H3K4me1 input
echo "H3K27ac_Tag_and_Signal"
echo "Input_Tag"
bam="MATERIAL/H3K4me1/merge/02.Align/H3K4me1_input_merge.bam"
input_tagdir="MATERIAL/H3K4me1/merge/03.TagDir/H3K4me1_input"
[[ -f "$bam" ]] && {
    mkdir -p "$input_tagdir"
    makeTagDirectory "$input_tagdir" "$bam"
} || echo "[WARNING] BAM not found for H3K4me1 input ($bam)"

# 6. H3K4me1
echo "H3K27ac_Tag_to_Signal"
bam="MATERIAL/H3K4me1/merge/02.Align/H3K4me1_merge.bam"
tagdir="MATERIAL/H3K4me1/merge/03.TagDir/H3K4me1"
[[ -f "$bam" ]] && {
    mkdir -p "$tagdir"
    makeTagDirectory "$tagdir" "$bam"

    outdir="MATERIAL/H3K4me1/merge/04.bigwig"
    mkdir -p "$outdir"
    makeUCSCfile "$tagdir" -bigWig "$CHROMSIZES" -c "$input_tagdir" -pseudo 1 -o "$outdir/H3K4me1.bw"
} || echo "[WARNING] BAM not found for H3K4me1 ($bam)"

echo ""
echo "All TagDirectories and bigWig files created."
