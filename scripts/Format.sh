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

# Process each annotation file ending with 1000bp_ann.txt or 5000bp_ann.txt
for file in "$BASEDIR"/*bp_ann.txt; do
    [ -f "$file" ] || continue
    base=$(basename "$file" .txt)
    outfile="$OUTPUT_DIR/${base}.tsv"

    echo "Formatting $file → $outfile"

    # Replace all spaces with tabs only between columns (preserving tab structure)
    # and ensure consistent formatting
    awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}' "$file" > "$outfile"
done

echo "All files reformatted and saved to: $OUTPUT_DIR"


# Input GTF file
GTF="/scratch/dr27977/OUTPUT/ANNOTATE/refann.gtf"
OUT="/scratch/dr27977/OUTPUT/ANNOTATE/Clean/refann_formatted.tsv"

# Output header
echo -e "Chr\tSource\tFeature\tStart\tEnd\tScore\tStrand\tFrame\tGeneID\tTranscriptID\tGeneName\tGeneBiotype\tTranscriptBiotype" > "$OUT"

# Extract and format the GTF
awk -F'\t' '
BEGIN { OFS="\t" }
{
    # Extract attributes
    match($9, /gene_id "([^"]+)"/, gid)
    match($9, /transcript_id "([^"]+)"/, tid)
    match($9, /gene_name "([^"]+)"/, gname)
    match($9, /gene_biotype "([^"]+)"/, gtype)
    match($9, /transcript_biotype "([^"]+)"/, ttype)

    print $1, $2, $3, $4, $5, $6, $7, $8,
          (gid[1] ? gid[1] : "NA"),
          (tid[1] ? tid[1] : "NA"),
          (gname[1] ? gname[1] : "NA"),
          (gtype[1] ? gtype[1] : "NA"),
          (ttype[1] ? ttype[1] : "NA")
}' "$GTF" >> "$OUT"

echo "GTF formatted to: $OUT"
