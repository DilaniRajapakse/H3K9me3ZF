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

# Input paths
BASEDIR="/scratch/dr27977/OUTPUT/ANNOTATE/Clean"
GTF="$BASEDIR/refann_formatted.tsv"
OUTDIR="$BASEDIR/PeakFeatureOverlap"
mkdir -p "$OUTDIR"

# Step 1: Build lookup from GTF: transcript_id + feature_type → chr, feature, start, end, gene_id, transcript_id, gene_symbol
awk -F'\t' 'BEGIN { OFS="\t" }
    FNR == 1 { next }  # Skip header
    {
        key = $10 "_" $3  # transcript_id + feature
        gtf_map[key] = $1 "\t" $3 "\t" $4 "\t" $5 "\t" $9 "\t" $10 "\t" $11
    }
    END {
        for (k in gtf_map) print k "\t" gtf_map[k]
    }
' "$GTF" > /tmp/gtf_lookup.tsv

# Step 2: Loop through ann.tsv files and match against GTF
for file in "$BASEDIR"/*.1000bp_ann.tsv "$BASEDIR"/*.5000bp_ann.tsv; do
    [[ -f "$file" ]] || continue
    base=$(basename "$file" .tsv)
    outfile="$OUTDIR/${base}_gtf_overlap.tsv"

    echo "Processing: $base → $outfile"

    awk -F'\t' -v src="$base" 'BEGIN {
        OFS = "\t"
        # Load GTF transcript-feature map
        while ((getline < "/tmp/gtf_lookup.tsv") > 0) {
            split($1, id_parts, "_")
            k = id_parts[1] "_" id_parts[2]
            lookup[k] = $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8
        }

        # Output header
        print "SourceFile", "PeakID", "Chr", "Feature", "Start", "End", "GeneID", "TranscriptID", "GeneSymbol"
    }

    FNR > 1 {
        transcript = $9
        feature = $8
        key = transcript "_" feature
        if (key in lookup) {
            print src, $1, lookup[key]
        }
    }' "$file" > "$outfile"
done

echo "All overlap tables (with headers and source filenames) saved to: $OUTDIR"