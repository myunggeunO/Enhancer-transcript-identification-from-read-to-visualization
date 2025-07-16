#!/bin/bash
set -e
THREADS=16

# 1. Merge ATAC replicates from *_noM_sorted.bam
echo "Merge ATAC replicates"
merge_dir="MATERIAL/ATAC/merge/02.Align"
mkdir -p "$merge_dir"

find MATERIAL/ATAC/rep{1,2}/02.Align -name "*_noM_sorted.bam" > ATAC_bams.txt

if [[ -s ATAC_bams.txt ]]; then
    sambamba merge -t $THREADS "$merge_dir/ATAC_merge.bam" $(cat ATAC_bams.txt)
else
    echo "[ERROR] No ATAC *_noM_sorted.bam files found"
fi

rm -f ATAC_bams.txt


# 2. Merge H3K4me1 ChIP replicates from *_rmdup.sorted.bam
echo "Merge H3K4me1 replicates"
merge_dir="MATERIAL/H3K4me1/merge/02.Align"
mkdir -p "$merge_dir"

find MATERIAL/H3K4me1/rep{1,2}/02.Align -name "*_rmdup.sorted.bam" > K4me1_bams.txt

if [[ -s K4me1_bams.txt ]]; then
    sambamba merge -t $THREADS "$merge_dir/H3K4me1_merge.bam" $(cat K4me1_bams.txt)
else
    echo "[ERROR] No H3K4me1 *_rmdup.sorted.bam files found"
fi

rm -f K4me1_bams.txt


# 3. Merge H3K4me1 Input replicates
echo "Merge H3K4me1 Input replicates"
find MATERIAL/H3K4me1/input/rep{1,2}/02.Align -name "*_rmdup.sorted.bam" > K4me1_input_bams.txt

if [[ -s K4me1_input_bams.txt ]]; then
    sambamba merge -t $THREADS "$merge_dir/H3K4me1_input_merge.bam" $(cat K4me1_input_bams.txt)
else
    echo "[ERROR] No H3K4me1 input *_rmdup.sorted.bam files found"
fi

rm -f K4me1_input_bams.txt

echo ""
echo "All merging complete"
