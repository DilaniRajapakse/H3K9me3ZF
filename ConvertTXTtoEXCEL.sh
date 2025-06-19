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
GTF="$BASEDIR/refann.gtf"
TE_BED="$BASEDIR/peaks/TEann_35_0.1filt.sorted.bed"
OUTDIR="$BASEDIR/peaks/tss_noTE_output"
mkdir -p "$OUTDIR/annotated_genes"

# Loop through 1kb annotated peak files
for annfile in "$BASEDIR"/peaks/ann/*1000bp_ann.txt; do
    base=$(basename "$annfile" .1000bp_ann.txt)

    # Convert to BED format
    awk -F'\t' '{print $2"\t"$3"\t"$4"\t"$1}' "$annfile" > "$OUTDIR/${base}_1kb_TSS_peaks.bed"

    # Remove TE-overlapping peaks
    bedtools intersect -v -a "$OUTDIR/${base}_1kb_TSS_peaks.bed" -b "$TE_BED" > "$OUTDIR/${base}_noTE.bed"

    # Annotate TE-free peaks
    annotatePeaks.pl "$OUTDIR/${base}_noTE.bed" danRer11 -gtf "$GTF" > "$OUTDIR/annotated_genes/${base}_noTE_annot.txt"

    # Extract Gene Symbol and Ensembl ID (unique)
    awk 'NR > 1 {print $2 "\t" $12}' "$OUTDIR/annotated_genes/${base}_noTE_annot.txt" | sort -u > "$OUTDIR/annotated_genes/${base}_noTE_genes.tsv"

    # Report how many genes had TE-free peaks near TSS
    count=$(wc -l < "$OUTDIR/annotated_genes/${base}_noTE_genes.tsv")
    echo "$base: $count genes with TE-free H3K9me3 peaks within 1kb of TSS"
done