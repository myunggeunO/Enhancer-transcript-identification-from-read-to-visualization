#!/bin/bash
set -e
GENOME="mm"

# 1. H3K27ac (merged tmerge, BAMPE, broad)
echo "Call peaks for H3K27ac"
mark="H3K27ac"
signal_bam="MATERIAL/$mark/tmerge/02.Align/${mark}_merge_filtered_rmdup.sorted.bam"
input_bam="MATERIAL/$mark/tmerge/02.Align/${mark}_input_merge_filtered_rmdup.sorted.bam"
outdir="MATERIAL/$mark/peak_calling"
mkdir -p "$outdir"

macs3 callpeak -t "$signal_bam" -c "$input_bam" -f BAMPE -g $GENOME --broad -n "$mark" --outdir "$outdir"

# 2. ATAC-seq (rep1 + rep2, BAMPE, narrow)
echo "Call peaks for ATAC"
mark="ATAC"
bam1="MATERIAL/$mark/rep1/02.Align/SRR23648611_noM_sorted.bam"
bam2="MATERIAL/$mark/rep2/02.Align/SRR23648610_noM_sorted.bam"
outdir="MATERIAL/$mark/peak_calling"
mkdir -p "$outdir"

macs3 callpeak -t "$bam1" "$bam2" -f BAMPE -g $GENOME --nomodel --shift -100 --extsize 200 -q 0.01 -n "$mark" --outdir "$outdir"

# 3. H3K4me1 (rep1 and rep2 separately, BAMPE, broad)
for rep in rep1 rep2; do
  echo "Call peaks for H3K4me1 $rep"
  mark="H3K4me1"
  bam=$(find MATERIAL/$mark/$rep/02.Align -name "*_rmdup.sorted.bam" | head -n 1)
  input=$(find MATERIAL/$mark/input/$rep/02.Align -name "*_rmdup.sorted.bam" | head -n 1)
  outdir="MATERIAL/$mark/peak_calling/$rep"
  mkdir -p "$outdir"

  if [[ -f "$bam" ]]; then
    macs3 callpeak -t "$bam" -c "$input" -f BAM -g $GENOME --broad -n "${mark}_${rep}" --outdir "$outdir"
  else
    echo "[ERROR] No BAM found for $mark $rep"
  fi
done

# 4. Identification of H3K4me1 overlapped broadpeak
mark="H3K4me1"
peak_dir="MATERIAL/$mark/peak_calling"
rep1_peak="$peak_dir/rep1/${mark}_rep1_peaks.broadPeak"
rep2_peak="$peak_dir/rep2/${mark}_rep2_peaks.broadPeak"
overlap_dir="$peak_dir/overlapped_peak"
output_peak="$overlap_dir/${mark}_overlapped_peaks.broadpeak"

mkdir -p "$overlap_dir"

if [[ -f "$rep1_peak" && -f "$rep2_peak" ]]; then
  echo "Find overlapping peaks between rep1 and rep2"
  bedtools intersect -a "$rep1_peak" -b "$rep2_peak" > "$output_peak"
  echo "Overlapped peaks saved to: $output_peak"
else
  echo "[ERROR] One or both peak files are missing:"
  [[ ! -f "$rep1_peak" ]] && echo " - Missing: $rep1_peak"
  [[ ! -f "$rep2_peak" ]] && echo " - Missing: $rep2_peak"
fi

echo ""
echo "All peak calling complete."
