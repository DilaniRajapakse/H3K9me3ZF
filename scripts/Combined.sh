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

# Inputs
ANN_DIR="/scratch/dr27977/OUTPUT/ANNOTATE"
GTF_TSV="$ANN_DIR/Clean/refann_formatted.tsv"
TE_BED="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/TEann_35_0.1filt.bed"
OUTDIR="$ANN_DIR/with_te"
mkdir -p "$OUTDIR"

echo "Step 1: Build gene body BED from GTF..."
awk -F'\t' 'FNR > 1 {
    chr = $1; start = $4 - 1; end = $5; txid = $10;
    print chr "\t" start "\t" end "\t" txid "\t.\t."
}' "$GTF_TSV" | sort -k1,1 -k2,2n > /tmp/gene_coords.bed

echo "Step 2: Intersect genes with TE BED..."
bedtools intersect -a /tmp/gene_coords.bed -b "$TE_BED" -wa -wb > /tmp/gene_te_overlap.tsv

echo "Step 3: Summarize TE overlaps..."
awk 'BEGIN{FS=OFS="\t"}
{
    gene = $4;
    gene_len = $3 - $2;
    overlap_len = ($7 > $2 ? ($7 < $3 ? $7 : $3) : $3) - ($6 > $2 ? $6 : $2);
    if (overlap_len > 0) {
        total_overlap[gene] += overlap_len;
        tes[gene] = (tes[gene] ? tes[gene] ";" $8 : $8);
        gene_length[gene] = gene_len;
    }
}
END {
    for (g in gene_length) {
        pct = (total_overlap[g] / gene_length[g]) * 100;
        bin = (pct == 0) ? "0%" :
              (pct <= 5) ? "≤5%" :
              (pct <= 10) ? "≤10%" :
              (pct <= 25) ? "≤25%" :
              (pct <= 50) ? "≤50%" :
              (pct <= 75) ? "≤75%" : "100%";
        print g, sprintf("%.2f", pct), bin, tes[g];
    }
}' /tmp/gene_te_overlap.tsv > /tmp/gene_te_summary.tsv

echo "Step 4: Annotate each *bp_ann.tsv file with TE info..."

for file in "$ANN_DIR"/*.1000bp_ann.txt "$ANN_DIR"/*.5000bp_ann.txt; do
    [[ -f "$file" ]] || continue
    base=$(basename "$file" .txt)
    outfile="$OUTDIR/${base}_withTE.tsv"

    echo "→ Annotating $base"

    awk 'BEGIN {
        FS = OFS = "\t";
        while ((getline < "/tmp/gene_te_summary.tsv") > 0) {
            tx = $1;
            te_pct[tx] = $2;
            te_bin[tx] = $3;
            te_names[tx] = $4;
        }
    }
    NR == 1 {
        print $0, "TE_Overlap_Percent", "TE_Overlap_Bin", "Overlapping_TE_Names";
        next;
    }
    {
        txid = $9;
        print $0,
              (txid in te_pct ? te_pct[txid] : "0.00"),
              (txid in te_bin ? te_bin[txid] : "0%"),
              (txid in te_names ? te_names[txid] : "None");
    }' "$file" > "$outfile"
done

echo "Done. All output written to: $OUTDIR"