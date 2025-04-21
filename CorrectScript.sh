#!/bin/bash
#SBATCH --job-name=NewH3K9me3	                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1	                                # Single task job
#SBATCH --cpus-per-task=24		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=120gb			                                # Total memory for job
#SBATCH --time=50:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/dr27977/log.%j		    # Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --mail-user=dr27977@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)

OUTDIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published" 
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

BASEDIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published"

#ml STAR
#for file in $OUTDIR/*_R*.fastq.gz;
#do
# if [[ $prefix ]]; then
#        base=$(basename ${first} _R1.fastq.gz)
#         sh /home/dr27977/H3K9me3ZF/PE_trim_and_star.sh -o $OUTDIR -n $base -m one $first $file
#         prefix=
#     else
#         first=$file
#         prefix=${file%%_*}
#     fi
# done

 ##aligning to ecoli genome
# curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > $OUTDIR/ecoli_refseq.fa
 # note here that STAR suggests SAindex = 10 but that makes the alignment FAIL, do 8 instead
# STAR --runThreadN 20 --genomeSAindexNbases 8 --runMode genomeGenerate --genomeDir $OUTDIR/ecoli_genome --genomeFastaFiles $OUTDIR/ecoli_refseq.fa

# for file in $OUTDIR/trimmed2/*_val_*.fq.gz;
# do
#   if [[ $prefix ]]; then
#       base=$(basename ${first} _R1_val_1.fq.gz)
#       STAR --runThreadN 20 --genomeDir $OUTDIR/ecoli_genome --outFileNamePrefix $OUTDIR/bams/"$base"_ecoli \
#         --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
#         --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
 
#       STAR --runThreadN 20 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/bams/"$base"_ecoli \
#        --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
#        --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000

#       STAR --runThreadN 20 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/bams/"$base"_ecoli \
#        --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
#        --outMultimapperOrder Random --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000

#       STAR --runThreadN 20 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/bams/"$base"_ecoli \
#        --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
#        --outSAMprimaryFlag AllBestScore --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
#         prefix=
#     else
#         first=$file
#         prefix=${file%%_*}
#     fi
# done

# rm $OUTDIR/bams/${base}*SJ.out.tab

# if [ -d "$OUTDIR/bams/logs" ]
# then
#     mv $OUTDIR/bams/*Log* $OUTDIR/bams/logs
# else
#   mkdir $OUTDIR/bams/logs
#   mv $OUTDIR/bams/*Log* $OUTDIR/bams/logs
# fi

#module load SAMtools

# for file in $OUTDIR/bams/${base}*ecoliAligned.sortedByCoord.out.bam
# do
#   base=$(basename ${file} ecoliAligned.sortedByCoord.out.bam)
#   samtools view -bq1 $file | samtools sort - > $OUTDIR/bams/${base}_ecoli_q1.bam
# done

# for infile in $OUTDIR/bams/*_ecoli_q1.bam
# do
#  base=$(basename "$infile" _ecoli_q1.bam)
#   echo "$base total aligned reads -" >> $OUTDIR/bams/bam_stats.txt
#   samtools view -@ 24 -F 0x4 $OUTDIR/bams/${base}ecoliAligned.sortedByCoord.out.bam | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
#   echo "  $base total aligned reads (unique mappers) -" >> $OUTDIR/bams/bam_stats.txt
#   samtools view -@ 24 -F 0x4 $OUTDIR/bams/${base}ecoliAligned.sortedByCoord.out.bam | grep "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
#   echo "  $base total aligned reads (multi mappers) -" >> $OUTDIR/bams/bam_stats.txt
#   samtools view -@ 24 -F 0x4 $OUTDIR/bams/${base}ecoliAligned.sortedByCoord.out.bam | grep -v "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
#   echo "$base q1 aligned reads -" >> $OUTDIR/bams/bam_stats.txt
#   samtools view -@ 24 -F 0x4 $infile | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
#   echo "  $base q1 aligned reads (unique mappers) -" >> $OUTDIR/bams/bam_stats.txt
#   samtools view -@ 24 -F 0x4 $infile | grep "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
#   echo "  $base q1 aligned reads (multi mappers) -" >> $OUTDIR/bams/bam_stats.txt
#   samtools view -@ 24 -F 0x4 $infile | grep -v "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
# done
#
#
# ###Remove PCR duplicates
module load picard
module load SAMtools

for infile in $OUTDIR/bams3/*_q1.bam
do
  base=$(basename ${infile} _q1.bam)

  #Add read groups
  java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=$OUTDIR/bams3/${base}_q1.bam \
    O=$OUTDIR/bams3/${base}_q1_rg.bam \
    RGID=1 \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=${base}

  # Run MarkDuplicates on the BAM with read groups
  java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I=$OUTDIR/bams3/${base}_q1_rg.bam \
    O=$OUTDIR/bams3/${base}_nodups.bam \
    M=$OUTDIR/bams3/${base}_dupmetrics.txt \
    REMOVE_DUPLICATES=true
done

# #merging IgG samples from all time points to create uniformity in peak calling later