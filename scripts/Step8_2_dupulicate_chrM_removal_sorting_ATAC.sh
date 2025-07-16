#!/bin/bash
set -e
THREADS=16

# 1. Run fixmate for ATAC-seq rep1 and rep2
echo "Run fixmate for ATAC-seq rep1 and rep2"

find MATERIAL/ATAC/rep{1,2}/02.Align -type f -name "*_filtered.bam" | while read bam; do
    base=$(basename "$bam" _filtered.bam)
    align_dir=$(dirname "$bam")
    tmp_dir="${align_dir}/tmp_${base}"
    mid_dir="${align_dir}/mid_${base}"
    mkdir -p "$tmp_dir" "$mid_dir"

    echo
    echo "Processing $base"

    # Name sort
    sambamba sort -t "$THREADS" -n -o "${tmp_dir}/${base}_namesorted.bam" "$bam"

    # Fixmate
    samtools fixmate -@ "$THREADS" -m \
        "${tmp_dir}/${base}_namesorted.bam" \
        "${tmp_dir}/${base}_fixmate.bam"

    # Coordinate sort
    sambamba sort -t "$THREADS" -o \
        "${tmp_dir}/${base}_sorted.bam" \
        "${tmp_dir}/${base}_fixmate.bam"

    # 2. Run markdup for ATAC-seq rep1 and rep2
    echo "Run markdup for ATAC-seq rep1 and rep2"
    dedup_bam="${align_dir}/${base}_final.bam"
    sambamba markdup -t "$THREADS" -r \
        "${tmp_dir}/${base}_sorted.bam" \
        "$dedup_bam"

    # 3. Remove chrM from deduplicated BAM
    echo "Remove chrM from deduplicated BAM"
    keep_chroms="${tmp_dir}/keep_chroms.txt"
    samtools idxstats "$dedup_bam" | cut -f1 | grep -v -w 'chrM' | grep -v -w '*' > "$keep_chroms"
    chrM_removed_bam="${align_dir}/${base}_noM.bam"
    samtools view -@ "$THREADS" -b "$dedup_bam" $(cat "$keep_chroms") > "$chrM_removed_bam"

    # 4. Run final coordinate sort
    echo "Run final coordinate sort"
    final_clean_bam="${chrM_removed_bam%.bam}_sorted.bam"
    sambamba sort -t "$THREADS" -o "$final_clean_bam" "$chrM_removed_bam"

    # 5. Move intermediate files
    echo "Move intermediate files"
    mv "${tmp_dir}"/* "$mid_dir/"
    rmdir "$tmp_dir"
    rm -f "$chrM_removed_bam" "$keep_chroms" "$dedup_bam" "${dedup_bam}.bai"

    echo "Done: $final_clean_bam"
done

echo ""
echo "All ATAC-seq replicates processing complete"
