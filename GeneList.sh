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

# Load required modules
module load BEDTools/2.31.0-GCC-12.3.0
module load Homer/5.1-foss-2023a-R-4.3.2

# Define paths
BASEDIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published"
IN_DIR="$BASEDIR/peaks"
TE_BED="$IN_DIR/TEann_35_0.1filt.sorted.bed"
GTF="$BASEDIR/refann.gtf"
OUT_DIR="$IN_DIR/TE_filtered"
mkdir -p "$OUT_DIR"

echo "Step 1: Removing peaks that overlap TE ≥50%"
for infile in "$IN_DIR"/*final.bed; do
  base=$(basename "$infile" .final.bed)
  bedtools intersect -a "$infile" -b "$TE_BED" -f 0.50 -v > "$OUT_DIR/${base}_TEann_final.bed"
done

echo "Step 2: Annotating filtered peaks with HOMER"
for infile in "$OUT_DIR"/*_TEann_final.bed; do
  base=$(basename "$infile" _TEann_final.bed)
  annotatePeaks.pl "$infile" danRer11 -gtf "$GTF" > "$OUT_DIR/${base}.TE_maskann.txt"
done

echo "Step 3: Filtering peaks within ±1kb of TSS"
for infile in "$OUT_DIR"/*.TE_maskann.txt; do
  base=$(basename "$infile" .TE_maskann.txt)
  awk -F'\t' 'sqrt($10*$10) <= 1000' "$infile" > "$OUT_DIR/${base}.1000bp_TEann.txt"
done

echo "Done. Final files saved to: $OUT_DIR"

echo "Step 4: Converting .1000bp_TEann.txt files to .csv format"
for txtfile in "$OUT_DIR"/*.1000bp_TEann.txt; do
  base=$(basename "$txtfile" .txt)
  csvfile="$OUT_DIR/${base}.csv"
  awk 'BEGIN {OFS=","} {print $0}' "$txtfile" > "$csvfile"
done
echo "CSV conversion complete."
