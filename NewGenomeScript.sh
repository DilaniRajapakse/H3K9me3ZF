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

#module load STAR
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

# for file in $OUTDIR/trimmed/*_val_*.fq.gz; do
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

# ###Remove PCR duplicates
module load picard 
module load SAMtools 

for infile in $OUTDIR/bams/*_q1.bam
do
  base=$(basename ${infile} _q1.bam)

  #Add read groups
  java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=$OUTDIR/bams/${base}_q1.bam \
    O=$OUTDIR/bams/${base}_q1_rg.bam \
    RGID=1 \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=${base}

  # Run MarkDuplicates on the BAM with read groups
  java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I=$OUTDIR/bams/${base}_q1_rg.bam \
    O=$OUTDIR/bams/${base}_nodups.bam \
    M=$OUTDIR/bams/${base}_dupmetrics.txt \
    REMOVE_DUPLICATES=true
done

# #merging IgG samples from all time points to create uniformity in peak calling later
##see if the other replicates exist
module load SAMtools 
samtools merge -f $BASEDIR/bams/IgG_nodups.bam \
  $BASEDIR/bams/IgG_2.5hpf_nodups.bam \
  $BASEDIR/bams/IgG_24hpf_nodups.bam \
  $BASEDIR/bams/IgG_4.5hpf_nodups.bam \
  $BASEDIR/bams/4hpf_IgG_1_nodups.bam \
  $BASEDIR/bams/4hpf_IgG_2_nodups.bam \
  $BASEDIR/bams/3hpf_IgG_1_nodups.bam \
  $BASEDIR/bams/3hpf_IgG_2_nodups.bam \
  $BASEDIR/bams/3.5hpf_IgG_1_nodups.bam \
  $BASEDIR/bams/2hpf_IgG_1_nodups.bam \

samtools merge -f $BASEDIR/bams/IgG_ecoli_nodups.bam \
 $BASEDIR/bams/IgG_2.5hpf__ecoli*nodups.bam \
 $BASEDIR/bams/IgG_24hpf__ecoli*nodups.bam \
 $BASEDIR/bams/IgG_4.5hpf__ecoli*nodups.bam \
 $BASEDIR/bams/4hpf_IgG_1__ecoli*nodups.bam \
 $BASEDIR/bams/4hpf_IgG_2__ecoli*nodups.bam \
 $BASEDIR/bams/3hpf_IgG_1__ecoli*nodups.bam \
 $BASEDIR/bams/3hpf_IgG_2__ecoli*nodups.bam \
 $BASEDIR/bams/3.5hpf_IgG_1__ecoli*nodups.bam \
 $BASEDIR/bams/2hpf_IgG_1__ecoli*nodups.bam \

##Now we need to extract all the aligned reads in preperation for spike in normalization
module load BEDTools 

for infile in $BASEDIR/bams/*nodups.bam
do
  base=$(basename ${infile} .bam)
  bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $BASEDIR/bams/$base.btb.bed
done

##spike in normalization
mkdir $BASEDIR/bdgrphs

for file in $BASEDIR/bams/*.btb.bed;
do
  if [[ $prefix ]]; then
        base=$(basename ${file} _nodups.btb.bed)
        sh /home/dr27977/H3K9me3ZF/DNA_spike.kd.sh $file $first \
        100000 bga $BASEDIR/genome/chrNameLength.txt 1 1000 $BASEDIR/bdgrphs/"$base".norm.bga
        prefix=
    else
        first=$file
        prefix=${file%%_*}
    fi
done

module load Homer 
mkdir $BASEDIR/peaks

#Convert normalized bga to BED format
for infile in $BASEDIR/bdgrphs/*.norm.bga; do
  base=$(basename "$infile" .norm.bga)
  awk '{print $1 "\t" $2 "\t" $3 "\t" "+" "\t" "+" "\t" "+"}' "$infile" > $BASEDIR/peaks/$base.bgato.bed
done

#Make tag directories for each BED
for infile in $BASEDIR/peaks/*.bgato.bed; do
  base=$(basename "$infile" .bgato.bed)
  makeTagDirectory $BASEDIR/peaks/$base.BtB.tagdir "$infile" -format bed
done

#Call peaks using merged IgG control
for infile in $BASEDIR/peaks/*K9*.BtB.tagdir; do
  base=$(basename "$infile" .BtB.tagdir)
  findPeaks "$infile" -style histone -minDist 1000 -gsize 1.5e9 -F 4 \
    -i $BASEDIR/peaks/IgG.BtB.tagdir \
    -o $BASEDIR/peaks/$base.txt
done

#Convert HOMER output to BED
for infile in $BASEDIR/peaks/*.txt; do
  base=$(basename "$infile" .txt)
  sed '/^#/d' "$infile" | \
    awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $8 "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' | \
    sed 's/\.000000//g' > $BASEDIR/peaks/$base.peaks.bed
done

module load ChIP-R 
chipr -i $BASEDIR/peaks/2hpf_K9_1.peaks.bed $BASEDIR/peaks/2hpf_K9_2.peaks.bed -m 2 -o $BASEDIR/peaks/2hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/K9abcam_2.5hpf_1.peaks.bed $BASEDIR/peaks/K9abcam_2.5hpf_2.peaks.bed $BASEDIR/peaks/K9abcam_2.5hpf_3.peaks.bed -m 2 -o $BASEDIR/peaks/2.5hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/3hpf_K9_1.peaks.bed $BASEDIR/peaks/3hpf_K9_2.peaks.bed $BASEDIR/peaks/3hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/3hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/3.5hpf_K9_1.peaks.bed $BASEDIR/peaks/3.5hpf_K9_2.peaks.bed $BASEDIR/peaks/3.5hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/3.5hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/4hpf_K9_1.peaks.bed $BASEDIR/peaks/4hpf_K9_2.peaks.bed $BASEDIR/peaks/4hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/4hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/K9abcam_4.5hpf_1.peaks.bed $BASEDIR/peaks/K9abcam_4.5hpf_2.peaks.bed $BASEDIR/peaks/K9abcam_4.5hpf_3.peaks.bed -m 2 -o $BASEDIR/peaks/4.5hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/K9abcam_24hpf_1.peaks.bed $BASEDIR/peaks/K9abcam_24hpf_2.peaks.bed $BASEDIR/peaks/K9abcam_24hpf_3.peaks.bed -m 2 -o $BASEDIR/peaks/24hpf_K9_repPeaks

###make a blacklist file
module load Homer 
makeTagDirectory $BASEDIR/peaks/IgG.BtB.tagdir $BASEDIR/peaks/IgG.bgato.bed -format bed
findPeaks $BASEDIR/peaks/IgG.BtB.tagdir -style factor -o $BASEDIR/peaks/IgG.txt
sed '/^#/d' $BASEDIR/peaks/IgG.txt | \
  awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" "1" "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' > $BASEDIR/peaks/blacklist.bed

ml BEDTools 

###intersect the peaks with the blacklist file to make sure we aren't looking at sticky regions before this step
for infile in $BASEDIR/peaks/*_repPeaks_all.bed
do
  base=$( basename ${infile} _repPeaks_all.bed )
  bedtools intersect -a $infile -b $BASEDIR/peaks/blacklist.bed -v > $BASEDIR/peaks/"$base"_final.bed
done



