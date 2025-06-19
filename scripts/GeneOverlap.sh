#!/bin/bash
#SBATCH --job-name=CreateTables	                            # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1	                                        # Single task job
#SBATCH --cpus-per-task=8		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=80gb			                                # Total memory for job
#SBATCH --time=72:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/dr27977/logs/output/log.%j.out
#SBATCH --error=/scratch/dr27977/logs/error/log.%j.err		# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --mail-user=dr27977@uga.edu                         # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=ALL                                     # Mail events (BEGIN, END, FAIL, ALL)

## intersect TE annotation file with K9 "gene" peak files to identify genic H3K9 enrichment
module load BEDTools/2.31.0-GCC-12.3.0

# Define paths
OUTDIR="/scratch/dr27977/OUTPUT/ANNOTATE"
PEAK_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks"
TE_BED="$PEAK_DIR/TEann_35_0.1filt.bed"
GTF="$OUTDIR/refann.gtf"

# Rename peak directories
FILTERED_DIR="$OUTDIR/filtered_peaks/TE_ann"
ANNOTATED_DIR="$OUTDIR/annotated_peaks/TE_ann"

mkdir -p "$FILTERED_DIR"
mkdir -p "$ANNOTATED_DIR"

echo "Step 1: Filter out peaks overlapping ≥50% with TE BED..."
for infile in "$PEAK_DIR"/*_final.bed; do
    base=$(basename "$infile" _final.bed)
    echo "→ Filtering $base"
    bedtools intersect -a "$infile" -b "$TE_BED" -f 0.50 -v > "$FILTERED_DIR/${base}_TEann_final.bed"
done

echo "Step 2: Annotate filtered peaks using GTF..."
for infile in "$FILTERED_DIR"/*_TEann_final.bed; do
    base=$(basename "$infile" _TEann_final.bed)
    echo "→ Annotating $base"
    annotatePeaks.pl "$infile" danRer11 -gtf "$GTF" > "$ANNOTATED_DIR/${base}.TE_maskann.txt"
done

echo "Step 3: Subset to ≤1kb and ≤5kb TSS windows..."
for infile in "$ANNOTATED_DIR"/*.TE_maskann.txt; do
    base=$(basename "$infile" .TE_maskann.txt)

    echo "→ Filtering 1kb: $base"
    awk -F'\t' 'sqrt($10*$10) <=1000' "$infile" > "$ANNOTATED_DIR/${base}.1000bp_TEann.txt"

    echo "→ Filtering 5kb: $base"
    awk -F'\t' 'sqrt($10*$10) <=5000' "$infile" > "$ANNOTATED_DIR/${base}.5000bp_TEann.txt"
done

echo "Step 4: Convert 1kb and 5kb annotations to BED..."
for file in "$ANNOTATED_DIR"/*.1000bp_TEann.txt "$ANNOTATED_DIR"/*.5000bp_TEann.txt; do
    [[ -f "$file" ]] || continue
    base=$(basename "$file" .txt)
    echo "→ Converting $base to BED"
    cut -f2-4 "$file" > "$ANNOTATED_DIR/${base}.bed"
done

echo "Step 5: Final TE exclusion confirmation for BED files..."
for bedfile in "$ANNOTATED_DIR"/*.bed; do
    base=$(basename "$bedfile" .bed)
    echo "→ Checking TE exclusion: $base"
    bedtools intersect -a "$bedfile" -b "$TE_BED" -v > "$ANNOTATED_DIR/${base}_final.bed"
done

echo "Done. Results in:"
echo "- $FILTERED_DIR"
echo "- $ANNOTATED_DIR"