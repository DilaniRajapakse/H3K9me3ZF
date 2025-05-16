#!/bin/bash
#SBATCH --job-name=NewGenomeScript                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1	                                # Single task job
#SBATCH --cpus-per-task=24		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=120gb			                                # Total memory for job
#SBATCH --time=72:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/dr27977/log.%j		    # Location of standard output and error log files 
#SBATCH --mail-user=dr27977@uga.edu                    # Where to send mail 
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

OUTDIR="/scratch/dr27977/H3K9me3_Zebrafish/NewGenome" 
#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

HOMEDIR="/home/dr27977/H3K9me3ZF"
if [ ! -d $HOMEDIR ]
then
    mkdir -p $HOMEDIR
fi

BASEDIR="/scratch/dr27977/H3K9me3_Zebrafish/NewGenome"

module load STAR
for file in $OUTDIR/*_R*.fastq.gz;
do
 if [[ $prefix ]]; then
        base=$(basename ${first} _R1.fastq.gz)
         sh /home/dr27977/H3K9me3ZF/PE_trim_and_star.sh -o $OUTDIR -n $base -m one $first $file
         prefix=
     else
         first=$file
         prefix=${file%%_*}
     fi
 done

 ##aligning to ecoli genome
 curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > $OUTDIR/ecoli_refseq.fa
 # note here that STAR suggests SAindex = 10 but that makes the alignment FAIL, do 8 instead
 STAR --runThreadN 20 --genomeSAindexNbases 8 --runMode genomeGenerate --genomeDir $OUTDIR/ecoli_genome --genomeFastaFiles $OUTDIR/ecoli_refseq.fa

 for file in $OUTDIR/trimmed/*_val_*.fq.gz;# do
   if [[ $prefix ]]; then
       base=$(basename ${first} _R1_val_1.fq.gz)
       STAR --runThreadN 20 --genomeDir $OUTDIR/ecoli_genome --outFileNamePrefix $OUTDIR/bams/"$base"_ecoli \
         --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
         --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
 
       STAR --runThreadN 20 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/bams/"$base"_ecoli \
        --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
        --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000

       STAR --runThreadN 20 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/bams/"$base"_ecoli \
        --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
        --outMultimapperOrder Random --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000

       STAR --runThreadN 20 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/bams/"$base"_ecoli \
        --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
        --outSAMprimaryFlag AllBestScore --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
         prefix=
     else
         first=$file
         prefix=${file%%_*}
     fi
 done

