#!/bin/bash
set -e
THREADS=5
ref_index="reference_index/mm10"

# 1. Run alignment for H3K27ac (paired-end)
echo "Run alignment for H3K27ac (paired-end)"
find MATERIAL/H3K27ac -type d -name "01.Clean" | while read clean_dir; do
    align_dir="$(dirname "$clean_dir")/02.Align"
    mkdir -p "$align_dir"

    echo "Processing $clean_dir"
    for val1 in "$clean_dir"/*_1_val_1.fq.gz; do
        base=$(basename "$val1" _1_val_1.fq.gz)
        val2="${clean_dir}/${base}_2_val_2.fq.gz"

        if [[ -f "$val1" && -f "$val2" ]]; then
            echo "Aligning: $base"
            bowtie2 -p "$THREADS" -x "$ref_index" \
                -1 "$val1" -2 "$val2" -S "$align_dir/${base}.sam"

            if [[ -f "$align_dir/${base}.sam" ]]; then
                samtools view -b -q 10 "$align_dir/${base}.sam" > "$align_dir/${                                                                                                                                     base}_filtered.bam"
                rm "$align_dir/${base}.sam"
            else
                echo "[ERROR] SAM not created: $base"
            fi
        fi
    done
done

# 2. Run alignment for H3K4me1 (single-end)
echo "Run alignment for H3K4me1 (single-end)"
find MATERIAL/H3K4me1 -type d -name "01.Clean" | while read clean_dir; do
    align_dir="$(dirname "$clean_dir")/02.Align"
    mkdir -p "$align_dir"

    echo "Processing $clean_dir"
    for fq in "$clean_dir"/*.fq.gz; do
        base=$(basename "$fq" _trimmed.fq.gz)

        echo "Aligning: $base"
        bowtie2 -p "$THREADS" -x "$ref_index" -U "$fq" -S "$align_dir/${base}.sa                                                                                                                                     m"

        if [[ -f "$align_dir/${base}.sam" ]]; then
            samtools view -b -q 10 "$align_dir/${base}.sam" > "$align_dir/${base                                                                                                                                     }_filtered.bam"
            rm "$align_dir/${base}.sam"
        else
            echo "[ERROR] SAM not created: $base"
        fi
    done
done

# 3. Run alignment for ATAC (paired-end)
echo "Run alignment for ATAC (paired-end)"
find MATERIAL/ATAC -type d -name "01.Clean" | while read clean_dir; do
    align_dir="$(dirname "$clean_dir")/02.Align"
    mkdir -p "$align_dir"

    echo "Processing $clean_dir"
    for val1 in "$clean_dir"/*_1_val_1.fq.gz; do
        base=$(basename "$val1" _1_val_1.fq.gz)
        val2="${clean_dir}/${base}_2_val_2.fq.gz"

        if [[ -f "$val1" && -f "$val2" ]]; then
            echo "Aligning: $base"
            bowtie2 -p "$THREADS" --very-sensitive -X 2000 -x "$ref_index" \
                -1 "$val1" -2 "$val2" -S "$align_dir/${base}.sam"

            if [[ -f "$align_dir/${base}.sam" ]]; then
                samtools view -b -q 30 "$align_dir/${base}.sam" > "$align_dir/${                                                                                                                                     base}_filtered.bam"
                rm "$align_dir/${base}.sam"
            else
                echo "[ERROR] SAM not created: $base"
            fi
        fi
    done
done

# 4. Run alignment for GRO (single-end)
echo "Run alignment for GRO (single-end)"
find MATERIAL/GRO -type d -name "01.Clean" | while read clean_dir; do
    align_dir="$(dirname "$clean_dir")/02.Align"
    mkdir -p "$align_dir"

    echo "Processing $clean_dir"
    for fq in "$clean_dir/"*_trimmed_polyA.fq.gz; do
        base=$(basename "$fq" _trimmed_polyA.fq.gz)

        echo "Aligning: $base"
        bowtie2 -p "$THREADS" --very-sensitive -x "$ref_index" -U "$fq" -S "$ali                                                                                                                                     gn_dir/${base}.sam"

        if [[ -f "$align_dir/${base}.sam" ]]; then
            samtools view -b -q 30 "$align_dir/${base}.sam" > "$align_dir/${base                                                                                                                                     }_filtered.bam"
            rm "$align_dir/${base}.sam"
        else
            echo "[ERROR] SAM not created: $base"
        fi
    done
done

echo ""
echo "All alignments complete"
