#!/bin/bash

##set expected directories for base script, input, and output
BASE_DIR="/scratch/dr27977"
INPUT_DIR="${BASE_DIR}/H3K9me3_Zebrafish/CUTnRUN_published/peaks"

##  /scratch/dr27977/OUTPUT/ANNOTATE

OUTPUT_DIR="${BASE_DIR}/OUTPUT/ANNOTATE"

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

##script start

##load modules
module load Homer/5.1-foss-2023a-R-4.3.2

##get reference annotation
curl -s ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz | gunzip -c > $OUTPUT_DIR/refann.gtf

for infile in ${INPUT_DIR}/*final.bed
  do
    base=$( basename ${infile} _final.bed)
    annotatePeaks.pl $infile danRer11 -gtf $OUTPUT_DIR/refann.gtf > $OUTPUT_DIR/${base}.maskann.txt
  done

##Filter peaks within 1kb of TSS
for infile in $OUTPUT_DIR/*maskann.txt
  do
    base=$(basename ${infile} .maskann.txt)
    awk -F'\t' 'sqrt($10*$10) <=1000' $infile > $OUTPUT_DIR/${base}.1000bp_ann.txt
  done

## Filter peaks greater than 1Kb of TSS
for infile in $OUTPUT_DIR/*maskann.txt
  do
    base=$(basename ${infile} .maskann.txt)
    awk -F'\t' 'sqrt($10*$10) >=1000' $infile | awk '{print $2 "\t" $3 "\t" $4 }' > $OUTPUT_DIR/${base}.MOREthan1000bp.bed
  done

##Filter peaks within 5kb of TSS
for infile in $OUTPUT_DIR/*maskann.txt
  do
    base=$(basename ${infile} .maskann.txt)
    awk -F'\t' 'sqrt($10*$10) <=5000' $infile > $OUTPUT_DIR/${base}.5000bp_ann.txt
  done

## Filter peaks greater than 5Kb of TSS
for infile in $OUTPUT_DIR/*maskann.txt
  do
    base=$(basename ${infile} .maskann.txt)
    awk -F'\t' 'sqrt($10*$10) >=5000' $infile | awk '{print $2 "\t" $3 "\t" $4 }' > $OUTPUT_DIR/${base}.MOREthan5000bp.bed
  done
