#!/bin/bash
#SBATCH --job-name=trim		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1		                            # Single task job
#SBATCH --cpus-per-task=24		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=120gb			                            # Total memory for job
#SBATCH --time=50:00:00  		                       # Time limit hrs:min:sec
#SBATCH --output=/scratch/dr27977/log.%j				# Location of standard output and error log files 
#SBATCH --mail-user=dr27977@uga.edu                    # Where to send mail 
#SBATCH --mail-type=BEGIN,END,FAIL,ALL                            # Mail events (BEGIN, END, FAIL, ALL)

OUTDIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_Abcam" 
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

BASEDIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_Abcam"
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
        #STAR --runThreadN 20 --genomeDir $BASEDIR/ecoli_genome --outFileNamePrefix $BASEDIR/bams2/"$base"_ecoli \
        #--readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
        #--outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
        #prefix=
    #else
        #first=$file
        #prefix=${file%%_*}
    #fi
#done

 
#for file in $OUTDIR/trimmed/*_val_*.fq.gz;
#do
#if [[ $prefix ]]; then
        #base=$(basename ${first} _R1_val_1.fq.gz)
        #STAR --runThreadN 20 --genomeDir $OUTDIR/ecoli_genome --outFileNamePrefix $OUTDIR/bams2/"$base"_ecoli \
        #--readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
        #--outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000

        #STAR --runThreadN 20 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/bams2/"$base"_ecoli \
        #--readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
        #--outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000

        #STAR --runThreadN 20 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/bams2/"$base"_ecoli \
        #--readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
        #--outMultimapperOrder Random --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000

        #STAR --runThreadN 20 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/bams2/"$base"_ecoli \
        #--readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
        #--outSAMprimaryFlag AllBestScore --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
        #prefix=
        #else
          #first=$file
          #prefix=${file%%_*}
       #fi
       #done
#Begin 4.9.25
#rm $OUTDIR/bams2/"$base"*SJ.out.tab

 #if [ -d "$OUTDIR/bams2/logs" ]
 #then
     #mv $OUTDIR/bams2/*Log* $OUTDIR/bams2/logs
 #else
   #mkdir $OUTDIR/bams2/logs
   #mv $OUTDIR/bams2/*Log* $OUTDIR/bams2/logs
 #fi

#module load SAMtools/1.18-GCC-12.3.0 

#for file in $OUTDIR/bams2/"$base"*ecoliAligned.sortedByCoord.out.bam
#do
  #base=$(basename ${file} ecoliAligned.sortedByCoord.out.bam)
  #samtools view -bq1 $file | samtools sort - > $OUTDIR/bams2/"$base"_ecoli_q1.bam
#done

#for infile in $OUTDIR/bams2/"$base"*_ecoli_q1.bam
#do
 #base=$(basename ${infile} _ecoli_q1.bam)
 #echo "$base total aligned reads -" >> $OUTDIR/bams2/bam_stats.txt
 #samtools view -@ 24 -F 0x4 $OUTDIR/bams2/"$base"ecoliAligned.sortedByCoord.out.bam | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams2/bam_stats.txt
 #echo "  $base total aligned reads (unique mappers) -" >> $OUTDIR/bams2/bam_stats.txt
 #samtools view -@ 24 -F 0x4 $OUTDIR/bams2/"$base"ecoliAligned.sortedByCoord.out.bam | grep "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams2/bam_stats.txt
 #echo "  $base total aligned reads (multi mappers) -" >> $OUTDIR/bams2/bam_stats.txt
 #samtools view -@ 24 -F 0x4 $OUTDIR/bams2/"$base"ecoliAligned.sortedByCoord.out.bam | grep -v "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams2/bam_stats.txt
 #echo "$base q1 aligned reads -" >> $OUTDIR/bams2/bam_stats.txt
 #samtools view -@ 24 -F 0x4 $infile | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams2/bam_stats.txt
 #echo "  $base q1 aligned reads (unique mappers) -" >> $OUTDIR/bams2/bam_stats.txt
 #samtools view -@ 24 -F 0x4 $infile | grep "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams2/bam_stats.txt
 #echo "  $base q1 aligned reads (multi mappers) -" >> $OUTDIR/bams2/bam_stats.txt
 #samtools view -@ 24 -F 0x4 $infile | grep -v "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams2/bam_stats.txt
 #done

#Begin 4.10.25
####Remove PCR duplicates
#ml picard/3.2.0-Java-17
module load SAMtools/1.18-GCC-12.3.0

#Check to see if read group information is present in a file. No read groups made
#samtools view -H $BASEDIR/bams/K9abcam_4.5hpf_3_ecoliAligned.sortedByCoord.out.bam | grep '@RG' >> $BASEDIR/bams/read_groups.txt

for infile in $OUTDIR/bams2/"$base"*_ecoliAligned.sortedByCoord.out.bam
do
  base=$(basename ${infile} _ecoliAligned.sortedByCoord.out.bam)
  samtools view -H $OUTDIR/bams2/"$base"_ecoliAligned.sortedByCoord.out.bam | grep '@RG' >> $OUTDIR/bams2/logs/read_groups.txt
done

#for infile in $BASEDIR/bams2/*q1.bam
#do
  #base=$(basename ${infile} _q1.bam)
  #java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $BASEDIR/bams2/"$base"_dupmetrics.txt -O $BASEDIR/bams2/"$base"_nodups.bam --REMOVE_DUPLICATES true
#done

 