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

BASE_DIR="/scratch/dr27977/OUTPUT/ANNOTATE"
OUTPUT_DIR="/scratch/dr27977/OUTPUT/ANNOTATE/Clean"
mkdir -p "$OUTPUT_DIR"

# Correct pattern based on actual filenames
for file in "$BASE_DIR"/*.1000bp_ann.txt "$BASE_DIR"/*.5000bp_ann.txt; do
    [ -f "$file" ] || continue
    base=$(basename "$file" .txt)
    outfile="$OUTPUT_DIR/${base}.tsv"

    echo "Formatting $file → $outfile"

    # Format and write only rows with 19+ fields
    awk 'BEGIN {OFS="\t"} NF >= 19 {
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,
              $10,$11,$12,$13,$14,$15,$16,$17,$18,$19
    }' "$file" > "$outfile"
done

echo "All files reformatted and saved to: $OUTPUT_DIR"