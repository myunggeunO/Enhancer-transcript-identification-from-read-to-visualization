#!/bin/bash
set -e

# 1. Set number of threads for fasterq-dump and pigz
THREADS=5

# 2. Create directory structure
mkdir -p MATERIAL/ATAC/{rep1/00.Rawdata,rep2/00.Rawdata}
mkdir -p MATERIAL/H3K27ac/{trep1/00.Rawdata,trep2/00.Rawdata,input/trep1/00.Rawdata,input/trep2/00.Rawdata}
mkdir -p MATERIAL/H3K4me1/{rep1/00.Rawdata,rep2/00.Rawdata,input/rep1/00.Rawdata,input/rep2/00.Rawdata}
mkdir -p MATERIAL/GRO/rep1/00.Rawdata

# 3. Define sample list
sample_list=$(cat <<EOF
ES_H3K27ac_trep1	SRR7895919
ES_H3K27ac_trep2	SRR7895920
ES_H3K27ac_input_trep1	SRR7895931
ES_H3K27ac_input_trep2	SRR7895932
ES_H3K4me1_rep1	SRR10049470
ES_H3K4me1_rep2	SRR10049471
ES_H3K4me1_input_rep1	SRR10049539
ES_H3K4me1_input_rep2	SRR10049540
ES_ATAC_rep1	SRR23648611
ES_ATAC_rep2	SRR23648610
ES_GRO_rep1	SRR5655667
EOF
)

# 4. Process each sample
echo "$sample_list" | while IFS=$'\t' read -r name srr; do
    echo "Processing $srr ($name)"

    # Parse components from name
    category=$(echo "$name" | cut -d'_' -f2)
    rep=$(echo "$name" | cut -d'_' -f3)
    last=$(echo "$name" | rev | cut -d'_' -f1 | rev)

    # Build rep_folder
    if [[ "$rep" == "input" ]]; then
        rep_folder="input/$last"
    else
        rep_folder="$rep"
    fi

    # Final path to 00.Rawdata
    target_dir="MATERIAL/$category/$rep_folder/00.Rawdata"
    echo "Target directory: $target_dir"

    if [[ ! -d "$target_dir" ]]; then
        echo "ERROR: Directory $target_dir does not exist"
        exit 1
    fi

    cd "$target_dir"

    echo "Downloading $srr with prefetch..."
    prefetch --output-directory . "$srr"

    sra_path="$PWD/$srr/$srr.sra"

    echo "Converting $srr to FASTQ..."
    if [[ "$category" == "ATAC" || "$category" == "H3K27ac" ]]; then
        fasterq-dump --split-files --threads $THREADS "$sra_path"
    else
        fasterq-dump --threads $THREADS "$sra_path"
    fi

    echo "Compressing FASTQ with pigz..."
    pigz -p $THREADS *.fastq

    echo "Cleaning up temporary files..."
    rm -rf "$srr"

    echo "$srr done."
    cd - > /dev/null
done

echo "All samples processed successfully."
