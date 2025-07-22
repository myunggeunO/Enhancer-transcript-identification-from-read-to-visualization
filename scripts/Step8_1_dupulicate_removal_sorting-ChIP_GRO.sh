#!/bin/bash
set -e
THREADS=25

# 1. Removing dupulicate and sorting for H3K27ac and corresponding inputs (for merged technical replicates)
echo "Remove dupulicate & sort for H3K27ac and inputs"
find MATERIAL/H3K27ac -path "*/tmerge/02.Align/*merge_filtered.bam" | while read bam; do
    base="${bam%.bam}"
    sambamba markdup -t $THREADS -r "$bam" "${base}_rmdup.bam"
    sambamba sort -t $THREADS -o "${base}_rmdup.sorted.bam" "${base}_rmdup.bam"
    rm -f "${base}_rmdup.bam" "${base}_rmdup.bam.bai"
done

# 2. Removing dupulicate and sorting for H3K4me1 and corresponding inputs
echo "Remove dupulicate & sort for H3K4me1 and inputs"
find MATERIAL/H3K4me1 -type f -name "*_filtered.bam" | while read bam; do
    base="${bam%.bam}"
    sambamba markdup -t $THREADS -r "$bam" "${base}_rmdup.bam"
    sambamba sort -t $THREADS -o "${base}_rmdup.sorted.bam" "${base}_rmdup.bam"
    rm -f "${base}_rmdup.bam" "${base}_rmdup.bam.bai"
done

# 3. Sort only for GRO files (no markdup)
echo "Sorting GRO-seq reads"
find MATERIAL/GRO -type f -name "*_filtered.bam" | while read bam; do
    base="${bam%.bam}"
    sambamba sort -t $THREADS -o "${base}.sorted.bam" "$bam"
done

echo ""
echo "ChIP, GRO rmdup, sorting processing complete"
