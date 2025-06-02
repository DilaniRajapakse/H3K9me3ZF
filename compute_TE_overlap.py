import pandas as pd
import csv
import pybedtools
import sys

annot_file = sys.argv[1]
trimmed_bed = sys.argv[2]
output_tsv = sys.argv[3]
timepoint = sys.argv[4]
window = sys.argv[5]

te_bed_path = "/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/TEann_35_0.1filt.bed"
symbol_tsv = "/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/H3K9me3_summary_tables/with_symbols/zebrafish_ensid_to_symbol.tsv"

# Load HOMER annotation
annot = pd.read_csv(annot_file, sep="\t", comment="#")
annot.columns = annot.columns.str.strip()

# Fix columns
if "PeakID" in annot.columns:
    annot = annot.rename(columns={"PeakID": "peak_id"})
elif annot.columns[0].startswith("PeakID"):
    annot = annot.rename(columns={annot.columns[0]: "peak_id"})

annot = annot.rename(columns={
    "Gene Name": "gene_symbol",
    "Nearest Ensembl": "gene_id",
    "Annotation": "annotation"
})

annot = annot.dropna(subset=["gene_id"])

# Classify peak type
annot["multiple_exons"] = annot["annotation"].str.contains("exon") & ~annot["annotation"].str.contains("exon 1")
annot["first_exon"] = annot["annotation"].str.contains("exon 1")
annot["multiple_introns"] = annot["annotation"].str.contains("intron")

# Peak lengths
peaks = pybedtools.BedTool(trimmed_bed)
peak_lengths = {i.name: int(i.end) - int(i.start) for i in peaks}

# TE overlap
te = pybedtools.BedTool(te_bed_path)
overlaps = {}
for i in peaks.intersect(te, wo=True):
    pid = i.name
    overlaps[pid] = overlaps.get(pid, 0) + int(i.fields[-1])

annot["TE_bp"] = annot["peak_id"].map(overlaps).fillna(0)
annot["Peak_bp"] = annot["peak_id"].map(peak_lengths).fillna(1)
annot["TE_pct"] = (annot["TE_bp"] / annot["Peak_bp"]) * 100

bins = [-0.1, 0, 10, 25, 50, 75, 100]
labels = ["0%", "<=10%", "<=25%", "<=50%", "<=75%", "100%"]
annot["TE_bin"] = pd.cut(annot["TE_pct"], bins=bins, labels=labels, include_lowest=True)

annot["timepoint"] = timepoint
annot["window"] = window

# Merge gene symbols
symbols = pd.read_csv(symbol_tsv, sep="\t", header=None, names=["gene_id", "gene_symbol_sym"])
annot = pd.merge(annot, symbols, on="gene_id", how="left")
annot["gene_symbol"] = annot["gene_symbol"].combine_first(annot["gene_symbol_sym"])

# Output summary
cols = ["gene_symbol", "gene_id", "multiple_exons", "multiple_introns", "first_exon", "TE_bin", "timepoint", "window"]
annot[cols].drop_duplicates().to_csv(output_tsv, sep="\t", index=False, quoting=csv.QUOTE_NONE)
