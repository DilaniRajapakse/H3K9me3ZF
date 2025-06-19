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

# Define input and output directories
INPUT_DIR="/scratch/dr27977/OUTPUT/ANNOTATE"
OUTPUT_DIR="$INPUT_DIR/excel_outputs"
mkdir -p "$OUTPUT_DIR"

# Load python if required by your system (uncomment if needed)
 module load pandas/1.0.5-foss-2022a-Python-3.10.4


# Find all matching txt files and convert to Excel
find "$INPUT_DIR" -type f \( -name "*K9.1000bp_ann.txt" -o -name "*K9.5000bp_ann.txt" \) | while read txtfile; do
    base=$(basename "$txtfile" .txt)
    outfile="$OUTPUT_DIR/${base}.xlsx"

    echo "Converting: $txtfile --> $outfile"

    python3 - <<EOF
import pandas as pd
df = pd.read_csv("$txtfile", sep='\t', header=None, dtype=str)
df.to_excel("$outfile", index=False, header=False)
EOF

done