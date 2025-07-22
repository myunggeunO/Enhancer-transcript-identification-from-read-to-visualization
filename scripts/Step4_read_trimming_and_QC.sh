#!/bin/bash
set -e

# Set number of threads
THREADS=5

# 1. Pre-trimming FastQC for all raw FASTQ files
echo "Performing initial FastQC on raw FASTQ files"

find MATERIAL -type f -path "*/00.Rawdata/*.fastq.gz" | while read fq; do
    echo "Running FastQC on $fq"
    fastqc -t $THREADS "$fq"
done

echo "Initial FastQC completed. Proceeding to trimming steps"

# 2. GRO-seq trimming
find MATERIAL/GRO -type d -name "00.Rawdata" | while read raw_dir; do
    clean_dir="$(dirname "$raw_dir")/01.Clean"
    mkdir -p "$clean_dir"

    for fq in "$raw_dir"/*.fastq.gz; do
        base=$(basename "$fq" .fastq.gz)
        echo "GRO: trimming $base"

        # 2.1. adapter trimming (Trim Galore)
        trim_galore -j "$THREADS" --nextseq 20 --length 20 --fastqc \
            -o "$clean_dir" "$fq"
            
        # 2.2. polyA trimming (Cutadapt)
        trimmed_fq="${clean_dir}/${base}_trimmed.fq.gz"
        cutadapt -j "$THREADS" -a "A{15}" -m 20 \
            -o "${clean_dir}/${base}_trimmed_polyA.fq.gz" "$trimmed_fq"

        # 2.3. Final QC
        fastqc -t $THREADS "${clean_dir}/${base}_trimmed_polyA.fq.gz"

        # 2.4. Remove intermediate file
        rm -f "$trimmed_fq"
    done
done

# 3. ATAC trimming
find MATERIAL/ATAC -type d -name "00.Rawdata" | while read raw_dir; do
    clean_dir="$(dirname "$raw_dir")/01.Clean"
    mkdir -p "$clean_dir"

    for r1 in "$raw_dir"/*_1.fastq.gz; do
        base=$(basename "$r1" _1.fastq.gz)
        r2="${raw_dir}/${base}_2.fastq.gz"
        echo "ATAC: trimming $base"
        trim_galore -j "$THREADS" --paired --nextera --fastqc -o "$clean_dir" "$r1" "$r2"
    done
done

# 4. H3K27ac trimming
find MATERIAL/H3K27ac -type d -name "00.Rawdata" | while read raw_dir; do
    clean_dir="$(dirname "$raw_dir")/01.Clean"
    mkdir -p "$clean_dir"

    for r1 in "$raw_dir"/*_1.fastq.gz; do
        base=$(basename "$r1" _1.fastq.gz)
        r2="${raw_dir}/${base}_2.fastq.gz"
        echo "H3K27ac: trimming $base"
        trim_galore -j "$THREADS" --paired --fastqc -o "$clean_dir" "$r1" "$r2"
    done
done

# 5. H3K4me1 trimming
find MATERIAL/H3K4me1 -type d -name "00.Rawdata" | while read raw_dir; do
    clean_dir="$(dirname "$raw_dir")/01.Clean"
    mkdir -p "$clean_dir"

    for fq in "$raw_dir"/*.fastq.gz; do
        base=$(basename "$fq" .fastq.gz)
        echo "H3K4me1: trimming $base"
        trim_galore -j "$THREADS" --fastqc -o "$clean_dir" "$fq"
    done
done

echo ""
echo "All trimming completed successfully."
