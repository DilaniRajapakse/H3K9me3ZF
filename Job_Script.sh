#!/bin/bash
#SBATCH --job-name=trim		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1		                            # Single task job
#SBATCH --cpus-per-task=24		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=120gb			                            # Total memory for job
#SBATCH --time=50:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/home/dr27977/H3K9me3ZF/log.%j				# Location of standard output and error log files 
#SBATCH --mail-user=dr27977@uga.edu                    # Where to send mail 
#SBATCH --mail-type=BEGIN,END,FAIL,ALL                            # Mail events (BEGIN, END, FAIL, ALL)

OUTDIR="/home/dr27977/H3K9me3_Zebrafish/CUTnRUN_Abcam" 
#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

HOMEDIR="/home/dr27977/H3K9me3ZF"
if [ ! -d $HOMEDIR ]
then
    mkdir -p $HOMEDIR
fi

BASEDIR="/home/dr27977/H3K9me3_Zebrafish/CUTnRUN_Abcam"
ml STAR

#for file in $BASEDIR/*_R*.fastq.gz;
#do
#if [[ $prefix ]]; then
        #base=$(basename ${first} _R1.fastq.gz)
        #sh $HOMEDIR/PE_trim_and_star.sh -o $BASEDIR -n $base -m one $first $file
        #prefix=
    #else
        #first=$file
        #prefix=${file%%_*}
    #fi
#done

#Files are being made, needed to update all the modules to newer versions. However, the job only ended because time ran out to run it. Increasing time 


###aligning to ecoli genome
#curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > $BASEDIR/ecoli_refseq.fa
###note here that STAR suggests SAindex = 10 but that makes the alignment FAIL, do 8 instead
#STAR --runThreadN 20 --genomeSAindexNbases 8 --runMode genomeGenerate --genomeDir $BASEDIR/ecoli_genome --genomeFastaFiles $BASEDIR/ecoli_refseq.fa

#for file in $BASEDIR/trimmed/*_val_*.fq.gz;
#do
  #if [[ $prefix ]]; then
        #base=$(basename ${first} _R1_val_1.fq.gz)
        #STAR --runThreadN 20 --genomeDir $BASEDIR/ecoli_genome --outFileNamePrefix $BASEDIR/bams/"$base"_ecoli \
        #--readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
        #--outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
        #prefix=
    #else
        #first=$file
        #prefix=${file%%_*}
    #fi
#done

module load SAMtools/1.18-GCC-12.3.0 

for file in $OUTDIR/bams/"$base"*ecoliAligned.sortedByCoord.out.bam
do
  base=$(basename ${file} ecoliAligned.sortedByCoord.out.bam)
  samtools view -bq1 $file | samtools sort - > $OUTDIR/bams/"$base"_ecoli_q1.bam
done

for infile in $OUTDIR/bams/"$base"*_ecoli_q1.bam
do
 base=$(basename ${infile} _ecoli_q1.bam)
 echo "$base total aligned reads -" >> $OUTDIR/bams/bam_stats.txt
 samtools view -@ 24 -F 0x4 $OUTDIR/bams/"$base"ecoliAligned.sortedByCoord.out.bam | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
 echo "  $base total aligned reads (unique mappers) -" >> $OUTDIR/bams/bam_stats.txt
 samtools view -@ 24 -F 0x4 $OUTDIR/bams/"$base"ecoliAligned.sortedByCoord.out.bam | grep "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
 echo "  $base total aligned reads (multi mappers) -" >> $OUTDIR/bams/bam_stats.txt
 samtools view -@ 24 -F 0x4 $OUTDIR/bams/"$base"ecoliAligned.sortedByCoord.out.bam | grep -v "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
 echo "$base q1 aligned reads -" >> $OUTDIR/bams/bam_stats.txt
 samtools view -@ 24 -F 0x4 $infile | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
 echo "  $base q1 aligned reads (unique mappers) -" >> $OUTDIR/bams/bam_stats.txt
 samtools view -@ 24 -F 0x4 $infile | grep "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
 echo "  $base q1 aligned reads (multi mappers) -" >> $OUTDIR/bams/bam_stats.txt
 samtools view -@ 24 -F 0x4 $infile | grep -v "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
 done

####Remove PCR duplicates
#ml picard/3.2.0-Java-17
#module load SAMtools/1.18-GCC-12.3.0

#for infile in $BASEDIR/bams/*q1.bam
#do
  #base=$(basename ${infile} _q1.bam)
  #java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $BASEDIR/bams/"$base"_dupmetrics.txt -O $BASEDIR/bams/"$base"_nodups.bam --REMOVE_DUPLICATES true
#done

 