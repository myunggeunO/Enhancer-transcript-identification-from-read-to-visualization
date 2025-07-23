#!/bin/bash
set -e
THREADS=5

# 1. Sort all bam of technical replicates
find MATERIAL/H3K27ac -type d -name "02.Align" | while read align_dir; do
    for bam in "$align_dir"/*_filtered.bam; do
        if [[ -f "$bam" ]]; then
            sorted_bam="${bam%.bam}.sorted.bam"
            sambamba sort -t $THREADS -o "$sorted_bam" "$bam"
        fi
    done
done

# 2. Create merge output directories
mkdir -p MATERIAL/H3K27ac/tmerge/02.Align

# 3. Merge files for ChIP techinical replicates (in this scripts, 27ac)
sambamba merge -t $THREADS MATERIAL/H3K27ac/tmerge/02.Align/H3K27ac_merge_filtered.bam \
    MATERIAL/H3K27ac/trep1/02.Align/*_filtered.sorted.bam \
    MATERIAL/H3K27ac/trep2/02.Align/*_filtered.sorted.bam

# 4. Merge files for Input replicates
sambamba merge -t $THREADS MATERIAL/H3K27ac/tmerge/02.Align/H3K27ac_input_merge_filtered.bam \
    MATERIAL/H3K27ac/input/trep1/02.Align/*_filtered.sorted.bam \
    MATERIAL/H3K27ac/input/trep2/02.Align/*_filtered.sorted.bam

echo ""
echo "Merging of technical replicates complete"
