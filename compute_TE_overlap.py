import pandas as pd
import csv
import pybedtools

# Load annotation
annot = pd.read_csv("2.5hpf_5kb_annot.txt", sep="\t")

# Clean column names
annot.columns = annot.columns.str.strip()

# Debugging print statements
print("Annotation Columns:")
print(annot.columns.tolist())
print("Preview of first few rows:")
print(annot.head(5))
print("Unique Annotations:")
print(annot['annotation'].dropna().unique())

print(f"Non-null gene_id count: {annot['gene_id'].notnull().sum()}")
print(f"Non-null peak_id count: {annot['peak_id'].notnull().sum()}")

# Rename only if the columns exist
if "PeakID" in annot.columns:
    annot = annot.rename(columns={"PeakID": "peak_id"})
elif annot.columns[0].startswith("PeakID"):
    annot = annot.rename(columns={annot.columns[0]: "peak_id"})

annot = annot.rename(columns={
    "Gene Name": "gene_symbol",
    "Nearest Ensembl": "gene_id",
    "Annotation": "annotation"
})

# Filter missing gene_id
annot = annot.dropna(subset=["gene_id"])

# Classify peak location
annot["multiple_exons"] = annot["annotation"].str.contains("exon") & ~annot["annotation"].str.contains("exon 1")
annot["first_exon"] = annot["annotation"].str.contains("exon 1")
annot["multiple_introns"] = annot["annotation"].str.contains("intron")

# Calculate peak lengths
peaks_bed = pybedtools.BedTool("2.5hpf_5kb_trimmed.bed")
peak_lengths = {}
for i in peaks_bed:
    peak_lengths[i.name] = int(i.end) - int(i.start)

# TE overlap per peak
te_bed = pybedtools.BedTool("/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/TEann_35_0.1filt.bed")
overlaps = {}
for i in peaks_bed.intersect(te_bed, wo=True):
    pid = i.name
    overlaps[pid] = overlaps.get(pid, 0) + int(i.fields[-1])

# Assign TE overlap % and bin
annot["TE_bp"] = annot["peak_id"].map(overlaps).fillna(0)
annot["Peak_bp"] = annot["peak_id"].map(peak_lengths).fillna(1)
annot["TE_pct"] = (annot["TE_bp"] / annot["Peak_bp"]) * 100

bins = [-0.1, 0, 10, 25, 50, 75, 100]
labels = ["0%", "<=10%", "<=25%", "<=50%", "<=75%", "100%"]
annot["TE_bin"] = pd.cut(annot["TE_pct"], bins=bins, labels=labels, include_lowest=True)

# Add timepoint and window manually for now
annot["timepoint"] = "2.5hpf"
annot["window"] = "5kb"

# Merge in gene symbols
symbols = pd.read_csv("/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/H3K9me3_summary_tables/with_symbols/zebrafish_ensid_to_symbol.tsv", sep="\t", header=None, names=["gene_id", "gene_symbol"])
annot = pd.merge(annot, symbols, on="gene_id", how="left", suffixes=("", "_sym"))

# Output columns
cols = ["gene_symbol", "gene_id", "multiple_exons", "multiple_introns", "first_exon", "TE_bin", "timepoint", "window"]
out = annot[cols].drop_duplicates()

out.to_csv("2.5hpf_K9_TSS5000bp_TE_table_with_symbols_summary.tsv", sep="\t", index=False, quoting=csv.QUOTE_NONE)
