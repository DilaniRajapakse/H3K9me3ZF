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

module load BEDTools
module load Homer/5.1-foss-2023a-R-4.3.2

# Define paths
BASEDIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published"
PEAKDIR="$BASEDIR/peaks"
OUTDIR="$BASEDIR/peaks/ann_TEfiltered"
TE_BED="$BASEDIR/peaks/TEann_35_0.1filt.bed"
GTF="$BASEDIR/refann.gtf"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Loop through each *_final.bed file
for bedfile in "$PEAKDIR"/*_final.bed; do
  base=$(basename "$bedfile" _final.bed)

  echo "Processing $base..."

  # Step 1: Remove peaks overlapping TEs (≥50% overlap)
  bedtools intersect -a "$bedfile" -b "$TE_BED" -f 0.50 -v > "$OUTDIR/${base}_TEfree.bed"

  # Step 2: Annotate TE-free peaks using HOMER and GTF
  annotatePeaks.pl "$OUTDIR/${base}_TEfree.bed" danRer11 -gtf "$GTF" > "$OUTDIR/${base}_TEfree.maskann.txt"

  # Step 3: Filter for ≤1000 bp peaks
  awk -F'\t' 'sqrt($10*$10) <= 1000' "$OUTDIR/${base}_TEfree.maskann.txt" > "$OUTDIR/${base}.1000bp_TEann.txt"

  # Step 4: Convert TXT to CSV
  awk 'BEGIN{OFS=","} {print $0}' "$OUTDIR/${base}.1000bp_TEann.txt" > "$OUTDIR/${base}.1000bp_TEann.csv"

done

echo "All files processed."