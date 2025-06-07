# Paths
BASE_DIR="/scratch/dr27977"
INPUT_DIR="${BASE_DIR}/OUTPUT/ANNOTATE/"
GTF="${INPUT_DIR}/refann.gtf"

##  /scratch/dr27977/OUTPUT/GENE

OUTPUT_DIR="${BASE_DIR}/OUTPUT/GENE"
TEMP_DIR="${BASE_DIR}/OUTPUT/GENE/TEMP"

##create directories if they do not currently exist
if [ ! -d $BASE_DIR ]
then
    mkdir -p $BASE_DIR
fi
cd $BASE_DIR

if [ ! -d $INPUT_DIR ]
then
    mkdir -p $INPUT_DIR
fi
cd $INPUT_DIR

if [ ! -d $OUTPUT_DIR ]
then
    mkdir -p $OUTPUT_DIR
fi
cd $OUTPUT_DIR

then
    mkdir -p $TEMP_DIR
fi
cd $TEMP_DIR

# Step 1: Extract ENSDART ➝ ENSDARG + Symbol mapping from GTF
awk -F'\t' '$3 == "transcript" {
    match($9, /transcript_id "([^"]+)"/, tid);
    match($9, /gene_id "([^"]+)"/, gid);
    match($9, /gene_name "([^"]+)"/, gname);
    if (tid[1] && gid[1] && gname[1]) print tid[1] "\t" gid[1] "\t" gname[1];
}' "$GTF" | sort -u > /${TEMP_DIR}/transcript_to_gene_symbol.tsv

# Step 2: Loop through each .maskann.txt file
for file in "$INPUT_DIR"/*.maskann.txt; do
    base=$(basename "$file" .maskann.txt)
    outfile="$OUTPUT_DIR/${base}_GeneTable.tsv"
    echo "Processing: $base"

    # Step 3: Extract required columns + match ENSDARG and symbol
    awk 'BEGIN {
        FS=OFS="\t";

        # Read ENSDART ➝ ENSDARG + GeneSymbol mapping
        while ((getline < "/tmp/transcript_to_gene_symbol.tsv") > 0) {
            txid2gid[$1] = $2;
            txid2sym[$1] = $3;
        }

        print "Chr", "Start", "End", "GeneID", "GeneSymbol";
    }
    FNR > 1 {
        tx = $14;
        gene_id = (tx in txid2gid ? txid2gid[tx] : "NA");
        gene_sym = (tx in txid2sym ? txid2sym[tx] : "NA");
        print $1, $2, $3, gene_id, gene_sym;
    }' "$file" > "$outfile"
done

echo "All done. Tables saved to: $OUTDIR"