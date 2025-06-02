import pandas as pd
import csv
import pybedtools

# Set your file paths and metadata here
annot_file = "PATH/TO/your_annot.txt"
trimmed_bed = "PATH/TO/your_trimmed.bed"
te_bed_path = "PATH/TO/TEann_35_0.1filt.bed"
symbol_tsv = "PATH/TO/zebrafish_ensid_to_symbol.tsv"
timepoint = "3hpf"
window = "5kb"
output = "PATH/TO/output_summary.tsv"

# Load HOMER annotation file
annot = pd.read_csv(annot_file, sep="\t", comment="#")
annot.columns = annot.columns.str.strip()

# Fix HOMER's inconsistent header naming
if "PeakID" in annot.columns:
    annot = annot.rename(columns={"PeakID": "peak_id"})
elif annot.columns[0].startswith("PeakID"):
    annot = annot.rename(columns={annot.columns[0]: "peak_id"})

annot = annot.rename(columns={
    "Gene Name": "gene_symbol",
    "Nearest Ensembl": "gene_id",
    "Annotation": "annotation"
})

# Drop entries without a gene ID
annot = annot.dropna(subset=["gene_id"])

# Classify peak overlap with gene features
annot["multiple_exons"] = annot["annotation"].str.contains("exon") & ~annot["annotation"].str.contains("exon 1")
annot["first_exon"] = annot["annotation"].str.contains("exon 1")
annot["multiple_introns"] = annot["annotation"].str.contains("intron")

# Compute peak lengths
peaks = pybedtools.BedTool(trimmed_bed)
peak_lengths = {i.name: int(i.end) - int(i.start) for i in peaks}

# Compute transposable element (TE) overlap in base pairs
te = pybedtools.BedTool(te_bed_path)
overlaps = {}
for i in peaks.intersect(te, wo=True):
    pid = i.name
    overlaps[pid] = overlaps.get(pid, 0) + int(i.fields[-1])

# Calculate % overlap with TE
annot["TE_bp"] = annot["peak_id"].map(overlaps).fillna(0)
annot["Peak_bp"] = annot["peak_id"].map(peak_lengths).fillna(1)
annot["TE_pct"] = (annot["TE_bp"] / annot["Peak_bp"]) * 100

# Bin genes by TE overlap
bins = [-0.1, 0, 10, 25, 50, 75, 100]
labels = ["0%", "<=10%", "<=25%", "<=50%", "<=75%", "100%"]
annot["TE_bin"] = pd.cut(annot["TE_pct"], bins=bins, labels=labels, include_lowest=True)

# Annotate metadata
annot["timepoint"] = timepoint
annot["window"] = window

# Add gene symbols from external table
symbols = pd.read_csv(symbol_tsv, sep="\t", header=None, names=["gene_id", "gene_symbol_sym"])
annot = pd.merge(annot, symbols, on="gene_id", how="left")
annot["gene_symbol"] = annot["gene_symbol"].combine_first(annot["gene_symbol_sym"])

# Final selected columns
cols = ["gene_symbol", "gene_id", "multiple_exons", "multiple_introns", "first_exon", "TE_bin", "timepoint", "window"]
annot[cols].drop_duplicates().to_csv(output, sep="\t", index=False, quoting=csv.QUOTE_NONE)
