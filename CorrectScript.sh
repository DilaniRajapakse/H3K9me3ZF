#!/bin/bash
#SBATCH --job-name=BinningGenes	                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1	                                # Single task job
#SBATCH --cpus-per-task=24		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=120gb			                                # Total memory for job
#SBATCH --time=72:00:00  		                            # Time limit hrs:min:sec
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

#ml STAR (SAMtools/1.18-GCC-12.3.0 )
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

#module load SAMtools (SAMtools/1.18-GCC-12.3.0) 

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
#module load picard (picard/3.2.0-Java-17 )
#module load SAMtools (SAMtools/1.18-GCC-12.3.0 )

#for infile in $OUTDIR/bams/*_q1.bam
#do
#  base=$(basename ${infile} _q1.bam)

  #Add read groups
#  java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#    I=$OUTDIR/bams/${base}_q1.bam \
#    O=$OUTDIR/bams/${base}_q1_rg.bam \
#    RGID=1 \
#    RGLB=lib1 \
#    RGPL=ILLUMINA \
#    RGPU=unit1 \
#    RGSM=${base}

  # Run MarkDuplicates on the BAM with read groups
#  java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
#    I=$OUTDIR/bams/${base}_q1_rg.bam \
#    O=$OUTDIR/bams/${base}_nodups.bam \
#    M=$OUTDIR/bams/${base}_dupmetrics.txt \
#    REMOVE_DUPLICATES=true
#done

# #merging IgG samples from all time points to create uniformity in peak calling later
#see if the other replicates exist
#module load SAMtools (SAMtools/1.18-GCC-12.3.0 )
#samtools merge -f $BASEDIR/bams/IgG_nodups.bam \
#  $BASEDIR/bams/IgG_2.5hpf_nodups.bam \
#  $BASEDIR/bams/IgG_24hpf_nodups.bam \
#  $BASEDIR/bams/IgG_4.5hpf_nodups.bam \
#  $BASEDIR/bams/4hpf_IgG_1_nodups.bam \
#  $BASEDIR/bams/4hpf_IgG_2_nodups.bam \
#  $BASEDIR/bams/3hpf_IgG_1_nodups.bam \
#  $BASEDIR/bams/3hpf_IgG_2_nodups.bam \
#  $BASEDIR/bams/3.5hpf_IgG_1_nodups.bam \
#  $BASEDIR/bams/2hpf_IgG_1_nodups.bam \

#samtools merge -f $BASEDIR/bams/IgG_ecoli_nodups.bam \
# $BASEDIR/bams/IgG_2.5hpf__ecoli*nodups.bam \
# $BASEDIR/bams/IgG_24hpf__ecoli*nodups.bam \
# $BASEDIR/bams/IgG_4.5hpf__ecoli*nodups.bam \
# $BASEDIR/bams/4hpf_IgG_1__ecoli*nodups.bam \
# $BASEDIR/bams/4hpf_IgG_2__ecoli*nodups.bam \
# $BASEDIR/bams/3hpf_IgG_1__ecoli*nodups.bam \
# $BASEDIR/bams/3hpf_IgG_2__ecoli*nodups.bam \
# $BASEDIR/bams/3.5hpf_IgG_1__ecoli*nodups.bam \
# $BASEDIR/bams/2hpf_IgG_1__ecoli*nodups.bam \

#Now we need to extract all the aligned reads in preperation for spike in normalization
#module load BEDTools (pybedtools/0.9.1-foss-2023a)

#for infile in $BASEDIR/bams/*nodups.bam
#do
#  base=$(basename ${infile} .bam)
#  bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $BASEDIR/bams/$base.btb.bed
#done

##spike in normalization
#mkdir $BASEDIR/bdgrphs

#for file in $BASEDIR/bams/*.btb.bed;
#do
#  if [[ $prefix ]]; then
#        base=$(basename ${file} _nodups.btb.bed)
#        sh /home/dr27977/H3K9me3ZF/DNA_spike.kd.sh $file $first \
#        100000 bga $BASEDIR/genome/chrNameLength.txt 1 1000 $BASEDIR/bdgrphs/"$base".norm.bga
#        prefix=
#    else
#        first=$file
#        prefix=${file%%_*}
#    fi
#done

###peak calling: I believe this is specifically for classifying peaks within 1kb of a genic TSS and the TEann portion is to intersect
### peaks with TE annotation to keep only TEs that overlap with at least 50% of the peak, to find TE-associated H3K9me3.
### TEcounts2.bed counts how many time each TE type is marked by H3K9me3 in that file


#module load Homer (Homer/5.1-foss-2023a-R-4.3.2)
#mkdir $BASEDIR/peaks

#Convert normalized bga to BED format
#for infile in $BASEDIR/bdgrphs/*.norm.bga; do
#  base=$(basename "$infile" .norm.bga)
#  awk '{print $1 "\t" $2 "\t" $3 "\t" "+" "\t" "+" "\t" "+"}' "$infile" > $BASEDIR/peaks/$base.bgato.bed
#done

#Make tag directories for each BED
#for infile in $BASEDIR/peaks/*.bgato.bed; do
#  base=$(basename "$infile" .bgato.bed)
#  makeTagDirectory $BASEDIR/peaks/$base.BtB.tagdir "$infile" -format bed
#done

#Call peaks using merged IgG control
#for infile in $BASEDIR/peaks/*K9*.BtB.tagdir; do
#  base=$(basename "$infile" .BtB.tagdir)
#  findPeaks "$infile" -style histone -minDist 1000 -gsize 1.5e9 -F 4 \
#    -i $BASEDIR/peaks/IgG.BtB.tagdir \
#    -o $BASEDIR/peaks/$base.txt
#done

#Convert HOMER output to BED
#for infile in $BASEDIR/peaks/*.txt; do
#  base=$(basename "$infile" .txt)
#  sed '/^#/d' "$infile" | \
#    awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $8 "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' | \
#    sed 's/\.000000//g' > $BASEDIR/peaks/$base.peaks.bed
#done

#module load ChIP-R (ChIP-R/1.1.0-foss-2022a-Python-3.10.4)
#chipr -i $BASEDIR/peaks/2hpf_K9_1.peaks.bed $BASEDIR/peaks/2hpf_K9_2.peaks.bed -m 2 -o $BASEDIR/peaks/2hpf_K9_repPeaks
#chipr -i $BASEDIR/peaks/K9abcam_2.5hpf_1.peaks.bed $BASEDIR/peaks/K9abcam_2.5hpf_2.peaks.bed $BASEDIR/peaks/K9abcam_2.5hpf_3.peaks.bed -m 2 -o $BASEDIR/peaks/2.5hpf_K9_repPeaks
#chipr -i $BASEDIR/peaks/3hpf_K9_1.peaks.bed $BASEDIR/peaks/3hpf_K9_2.peaks.bed $BASEDIR/peaks/3hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/3hpf_K9_repPeaks
#chipr -i $BASEDIR/peaks/3.5hpf_K9_1.peaks.bed $BASEDIR/peaks/3.5hpf_K9_2.peaks.bed $BASEDIR/peaks/3.5hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/3.5hpf_K9_repPeaks
#chipr -i $BASEDIR/peaks/4hpf_K9_1.peaks.bed $BASEDIR/peaks/4hpf_K9_2.peaks.bed $BASEDIR/peaks/4hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/4hpf_K9_repPeaks
#chipr -i $BASEDIR/peaks/K9abcam_4.5hpf_1.peaks.bed $BASEDIR/peaks/K9abcam_4.5hpf_2.peaks.bed $BASEDIR/peaks/K9abcam_4.5hpf_3.peaks.bed -m 2 -o $BASEDIR/peaks/4.5hpf_K9_repPeaks
#chipr -i $BASEDIR/peaks/K9abcam_24hpf_1.peaks.bed $BASEDIR/peaks/K9abcam_24hpf_2.peaks.bed $BASEDIR/peaks/K9abcam_24hpf_3.peaks.bed -m 2 -o $BASEDIR/peaks/24hpf_K9_repPeaks

###make a blacklist file
#module load Homer (Homer/5.1-foss-2023a-R-4.3.2)
#makeTagDirectory $BASEDIR/peaks/IgG.BtB.tagdir $BASEDIR/peaks/IgG.bgato.bed -format bed
#findPeaks $BASEDIR/peaks/IgG.BtB.tagdir -style factor -o $BASEDIR/peaks/IgG.txt
#sed '/^#/d' $BASEDIR/peaks/IgG.txt | \
#  awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" "1" "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' > $BASEDIR/peaks/blacklist.bed

#ml BEDTools (pybedtools/0.9.1-foss-2023a)

###intersect the peaks with the blacklist file to make sure we aren't looking at sticky regions before this step
#for infile in $BASEDIR/peaks/*_repPeaks_all.bed
#do
#  base=$( basename ${infile} _repPeaks_all.bed )
#  bedtools intersect -a $infile -b $BASEDIR/peaks/blacklist.bed -v > $BASEDIR/peaks/"$base"_final.bed
#done

###peak annotation#### 4.25.25 This section is to get a gene list with H3K9me3 enriched peaks. bga files before will show all H3K9me3 enrichment across the genome, we need maskann list to find specific genes.

##mask ann is all annotated H3K9me3 peaks in zebrafish genome
##TE ann is all transposoable elements in genome with H3K9me3

#5.7.25 Get gene list of H3K9me3 annotated peaks within 5kb of TSS and peaks outside the 5kb region
#module load Homer/5.1-foss-2023a-R-4.3.2
#module load BEDtools
#curl -s ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz | gunzip -c > $BASEDIR/refann.gtf
#mkdir $BASEDIR/peaks/ann

#for infile in $BASEDIR/peaks/*final.bed
#do
#  base=$( basename ${infile} final.bed)
#  annotatePeaks.pl $infile danRer11 -gtf $BASEDIR/refann.gtf > $BASEDIR/peaks/ann/$base.maskann.txt
#done

#for infile in $BASEDIR/peaks/ann/*maskann.txt
#do
#  base=$(basename ${infile} .maskann.txt)
#  awk -F'\t' 'sqrt($10*$10) <=5000' $infile > $BASEDIR/peaks/ann/$base.5000bp_ann.txt
#done

#for infile in $OUTDIR/peaks/ann/*maskann.txt
# do
#   base=$(basename ${infile} .maskann.txt)
#   awk -F'\t' 'sqrt($10*$10) >=5000' $infile | awk '{print $2 "\t" $3 "\t" $4 }' > $BASEDIR/peaks/ann/${base}.MOREthan5000bp.bed
# done
###  awk -F'\t' 'sqrt($10*$10) <=1000' $infile > $BASEDIR/peaks/ann/$base.1000bp_ann.txt. The 1000 is to get within 1kb of a gene
#Did not have annotated TE file uploaded



###4.23.25. Going forward, I believe I want to take my originally created .bga files and make bigwigs from those to make timepoint
###tracks in IGV? I think I need to merge timepoints together to create 1 bw per timepoint and then have 1 merged IgG track that I can use to compare to all time points

#4.24.25
###lets make these bedgraphs into bigwigs for data visualization
#module load ucsc/443
#mkdir $BASEDIR/bws

#for infile in $BASEDIR/bdgrphs/*norm.bga
#do
#  base=$(basename ${infile} .norm.bga)
#  bedSort $infile $infile
#  bedGraphToBigWig $infile $BASEDIR/genome/chrNameLength.txt $BASEDIR/bws/$base.bw
#done

###lets do some broad comparisons to see what our data looks like before moving on
#module load deepTools (deepTools/3.5.2-foss-2022a)


#multiBigwigSummary bins -b $BASEDIR/bws/*[1-3].bw $BASEDIR/bws/*IgG*.bw -o $BASEDIR/bwreps_summ.npz -p 24
#plotCorrelation -in $BASEDIR/bwreps_summ.npz -c spearman -p heatmap -o $BASEDIR/timecourse_bwreps_summ_heatmap.pdf
#plotPCA -in $BASEDIR/bwreps_summ.npz -o $BASEDIR/timecourse_bwreps_summ_PCA.pdf

## 5.7.25 making bws of the mean of the replicates for profile plots
#module load deepTools/3.5.2-foss-2022a
#mkdir $OUTDIR/bws/bwscomp/
#bigwigCompare \
#  -b1 $OUTDIR/bws/2hpf_K9_1.bw \
#  -b2 $OUTDIR/bws/2hpf_K9_2.bw \
#  --operation add \
#  --scaleFactors 0.5:0.5 \
#  -bs 10 \
#  -p 20 \
#  -o $OUTDIR/bws/bwscomp/2hpf_K9_AVG.bw

# Step 1: Combine replicate 1 and 2 (each scaled to 1/3)
#bigwigCompare \
#  -b1 $OUTDIR/bws/K9abcam_2.5hpf_1.bw \
#  -b2 $OUTDIR/bws/K9abcam_2.5hpf_2.bw \
#  --operation add \
#  --scaleFactors 0.333:0.333 \
#  -bs 10 \
#  -p 20 \
#  -o $OUTDIR/bws/bwscomp/K9abcam_2.5hpf_1plus2.bw

# Step 2: Add replicate 3 to the combined 1+2 (again scaled to 1/3)
#bigwigCompare \
#  -b1 $OUTDIR/bws/bwscomp/K9abcam_2.5hpf_1plus2.bw \
#  -b2 $OUTDIR/bws/K9abcam_2.5hpf_3.bw \
#  --operation add \
#  --scaleFactors 1:0.333 \
#  -bs 10 \
#  -p 20 \
#  -o $OUTDIR/bws/bwscomp/K9abcam_2.5hpf_AVG.bw

# Step 1: Combine replicates 1 and 2, each scaled to 1/3
#bigwigCompare \
#  -b1 $OUTDIR/bws/3hpf_K9_1.bw \
#  -b2 $OUTDIR/bws/3hpf_K9_2.bw \
#  --operation add \
#  --scaleFactors 0.333:0.333 \
#  -bs 10 \
#  -p 20 \
#  -o $OUTDIR/bws/bwscomp/3hpf_K9_1plus2.bw

# Step 2: Add replicate 3, also scaled to 1/3
#bigwigCompare \
#  -b1 $OUTDIR/bws/bwscomp/3hpf_K9_1plus2.bw \
#  -b2 $OUTDIR/bws/3hpf_K9_3.bw \
#  --operation add \
#  --scaleFactors 1:0.333 \
#  -bs 10 \
#  -p 20 \
#  -o $OUTDIR/bws/bwscomp/3hpf_K9_AVG.bw

# Step 1: Combine replicates 1 and 2 (each scaled to 1/3)
#bigwigCompare \
#  -b1 $OUTDIR/bws/3.5hpf_K9_1.bw \
#  -b2 $OUTDIR/bws/3.5hpf_K9_2.bw \
#  --operation add \
#  --scaleFactors 0.333:0.333 \
#  -bs 10 \
#  -p 20 \
#  -o $OUTDIR/bws/bwscomp/3.5hpf_K9_1plus2.bw

# Step 2: Add replicate 3 (also scaled to 1/3)
#bigwigCompare \
#  -b1 $OUTDIR/bws/bwscomp/3.5hpf_K9_1plus2.bw \
#  -b2 $OUTDIR/bws/3.5hpf_K9_3.bw \
#  --operation add \
#  --scaleFactors 1:0.333 \
#  -bs 10 \
#  -p 20 \
#  -o $OUTDIR/bws/bwscomp/3.5hpf_K9_AVG.bw

# Step 1: Combine replicates 1 and 2 (each scaled to 1/3)
#bigwigCompare \
#  -b1 $OUTDIR/bws/4hpf_K9_1.bw \
#  -b2 $OUTDIR/bws/4hpf_K9_2.bw \
#  --operation add \
#  --scaleFactors 0.333:0.333 \
#  -bs 10 \
#  -p 20 \
#  -o $OUTDIR/bws/bwscomp/4hpf_K9_1plus2.bw

# Step 2: Add replicate 3 (also scaled to 1/3)
#bigwigCompare \
#  -b1 $OUTDIR/bws/bwscomp/4hpf_K9_1plus2.bw \
#  -b2 $OUTDIR/bws/4hpf_K9_3.bw \
#  --operation add \
#  --scaleFactors 1:0.333 \
#  -bs 10 \
#  -p 20 \
#  -o $OUTDIR/bws/bwscomp/4hpf_K9_AVG.bw

# Step 1: Combine replicates 1 and 2 (each scaled to 1/3)
#bigwigCompare \
#  -b1 $OUTDIR/bws/K9abcam_4.5hpf_1.bw \
#  -b2 $OUTDIR/bws/K9abcam_4.5hpf_2.bw \
#  --operation add \
#  --scaleFactors 0.333:0.333 \
#  -bs 10 \
#  -p 20 \
#  -o $OUTDIR/bws/bwscomp/K9abcam_4.5hpf_1plus2.bw

# Step 2: Add replicate 3 (also scaled to 1/3)
#bigwigCompare \
#  -b1 $OUTDIR/bws/bwscomp/K9abcam_4.5hpf_1plus2.bw \
#  -b2 $OUTDIR/bws/K9abcam_4.5hpf_3.bw \
#  --operation add \
#  --scaleFactors 1:0.333 \
#  -bs 10 \
#  -p 20 \
#  -o $OUTDIR/bws/bwscomp/K9abcam_4.5hpf_AVG.bw

# Step 1: Combine replicates 1 and 2 (each scaled to 1/3)
#bigwigCompare \
#  -b1 $OUTDIR/bws/K9abcam_24hpf_1.bw \
#  -b2 $OUTDIR/bws/K9abcam_24hpf_2.bw \
#  --operation add \
#  --scaleFactors 0.333:0.333 \
#  -bs 10 \
#  -p 20 \
#  -o $OUTDIR/bws/bwscomp/K9abcam_24hpf_1plus2.bw

# Step 2: Add replicate 3 (also scaled to 1/3)
#bigwigCompare \
#  -b1 $OUTDIR/bws/bwscomp/K9abcam_24hpf_1plus2.bw \
#  -b2 $OUTDIR/bws/K9abcam_24hpf_3.bw \
#  --operation add \
#  --scaleFactors 1:0.333 \
#  -bs 10 \
#  -p 20 \
#  -o $OUTDIR/bws/bwscomp/K9abcam_24hpf_AVG.bw

#5.7.25- 5.8.25
module load BEDTools/2.31.0-GCC-12.3.0
module load Homer/5.1-foss-2023a-R-4.3.2 

PEAKS=$OUTDIR/peaks
TE_ANN=/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/TEann_35_0.1filt.bed 
REF_ANN=$OUTDIR/refann.gtf

BIN_DIR=$OUTDIR/TE_overlap_bins
mkdir -p $BIN_DIR

for infile in $PEAKS/*final.bed
do
  base=$(basename $infile _final.bed)
  echo "Processing $base..."

  # Temp file to hold peaks not yet assigned to a bin
  remaining_peaks=$BIN_DIR/${base}_remaining.bed
  cp $infile $remaining_peaks

  # --- Bin 1: 0% overlap ---
  mkdir -p $BIN_DIR/${base}_bin_000
  bedtools intersect -a $remaining_peaks -b $TE_ANN -v > $BIN_DIR/${base}_bin_000/${base}_TEbin_000.bed
  bedtools intersect -a $remaining_peaks -b $TE_ANN -u > $BIN_DIR/tmp_used.bed
  bedtools intersect -a $remaining_peaks -b $BIN_DIR/tmp_used.bed -v > $BIN_DIR/tmp_next.bed
  mv $BIN_DIR/tmp_next.bed $remaining_peaks

  # --- Bin 2: >0% and ≤10% ---
  mkdir -p $BIN_DIR/${base}_bin_010
  bedtools intersect -a $remaining_peaks -b $TE_ANN -f 0.0001 -u > $BIN_DIR/tmp_overlap_gt0.bed
  bedtools intersect -a $remaining_peaks -b $TE_ANN -f 0.10 -v > $BIN_DIR/tmp_overlap_le10.bed
  bedtools intersect -a $BIN_DIR/tmp_overlap_gt0.bed -b $BIN_DIR/tmp_overlap_le10.bed > $BIN_DIR/${base}_bin_010/${base}_TEbin_010.bed
  bedtools intersect -a $remaining_peaks -b $BIN_DIR/${base}_bin_010/${base}_TEbin_010.bed -v > $BIN_DIR/tmp_next.bed
  mv $BIN_DIR/tmp_next.bed $remaining_peaks
  rm $BIN_DIR/tmp_overlap_gt0.bed $BIN_DIR/tmp_overlap_le10.bed

  # --- Bin 3: >10% and ≤25% ---
  mkdir -p $BIN_DIR/${base}_bin_025
  bedtools intersect -a $remaining_peaks -b $TE_ANN -f 0.10 -u > $BIN_DIR/tmp_gt10.bed
  bedtools intersect -a $remaining_peaks -b $TE_ANN -f 0.25 -v > $BIN_DIR/tmp_le25.bed
  bedtools intersect -a $BIN_DIR/tmp_gt10.bed -b $BIN_DIR/tmp_le25.bed > $BIN_DIR/${base}_bin_025/${base}_TEbin_025.bed
  bedtools intersect -a $remaining_peaks -b $BIN_DIR/${base}_bin_025/${base}_TEbin_025.bed -v > $BIN_DIR/tmp_next.bed
  mv $BIN_DIR/tmp_next.bed $remaining_peaks
  rm $BIN_DIR/tmp_gt10.bed $BIN_DIR/tmp_le25.bed

  # --- Bin 4: >25% and ≤50% ---
  mkdir -p $BIN_DIR/${base}_bin_050
  bedtools intersect -a $remaining_peaks -b $TE_ANN -f 0.25 -u > $BIN_DIR/tmp_gt25.bed
  bedtools intersect -a $remaining_peaks -b $TE_ANN -f 0.50 -v > $BIN_DIR/tmp_le50.bed
  bedtools intersect -a $BIN_DIR/tmp_gt25.bed -b $BIN_DIR/tmp_le50.bed > $BIN_DIR/${base}_bin_050/${base}_TEbin_050.bed
  bedtools intersect -a $remaining_peaks -b $BIN_DIR/${base}_bin_050/${base}_TEbin_050.bed -v > $BIN_DIR/tmp_next.bed
  mv $BIN_DIR/tmp_next.bed $remaining_peaks
  rm $BIN_DIR/tmp_gt25.bed $BIN_DIR/tmp_le50.bed

  # --- Bin 5: >50% and ≤75% ---
  mkdir -p $BIN_DIR/${base}_bin_075
  bedtools intersect -a $remaining_peaks -b $TE_ANN -f 0.50 -u > $BIN_DIR/tmp_gt50.bed
  bedtools intersect -a $remaining_peaks -b $TE_ANN -f 0.75 -v > $BIN_DIR/tmp_le75.bed
  bedtools intersect -a $BIN_DIR/tmp_gt50.bed -b $BIN_DIR/tmp_le75.bed > $BIN_DIR/${base}_bin_075/${base}_TEbin_075.bed
  bedtools intersect -a $remaining_peaks -b $BIN_DIR/${base}_bin_075/${base}_TEbin_075.bed -v > $BIN_DIR/tmp_next.bed
  mv $BIN_DIR/tmp_next.bed $remaining_peaks
  rm $BIN_DIR/tmp_gt50.bed $BIN_DIR/tmp_le75.bed

  # --- Bin 6: >75% ---
  mkdir -p $BIN_DIR/${base}_bin_100
  mv $remaining_peaks $BIN_DIR/${base}_bin_100/${base}_TEbin_100.bed

  # Annotate each bin
  for binfile in $BIN_DIR/${base}_bin_*/${base}_TEbin_*.bed
  do
    bindir=$(dirname "$binfile")
    binbase=$(basename $binfile .bed)
    annotatePeaks.pl $binfile danRer11 -gtf $GTF > $bindir/${binbase}.ann.txt
    awk -F'\t' 'sqrt($10*$10) <=5000' $bindir/${binbase}.ann.txt > $bindir/${binbase}.within5kb.txt
    cut -f2 $bindir/${binbase}.within5kb.txt | tail -n +2 | sort | uniq > $bindir/${binbase}_genes.txt
  done

done


#4.23.25 Original: all of this is TE specific, I don't need those
#TEFILE="$BASEDIR/peaks/TEann_35_0.1filt.bed"
#ANNDIR="$BASEDIR/peaks/ann"
#ANNOUTDIR="$BASEDIR/peaks/ann4"

#for infile in "$ANNDIR"/*maskann.txt; do
#  base=$(basename "$infile" .maskann.txt)
#  awk -F'\t' 'sqrt($10*$10) >= 1000' "$infile" | awk '{print $2 "\t" $3 "\t" $4 }' > "$ANNDIR/$base.MOREthan1000bp.bed"
#done

###Filtering for broad peaks (greater than or = to 1kb) using sqrt($10*$10) >= 1000 (1000bp peaks from HOMER’s annotation) So if I'm interested in 5kb region I have to change this
#mkdir -p "$ANNOUTDIR"

#for infile in "$ANNDIR"/*MOREthan1000bp.bed; do
#  base=$(basename "$infile" .MOREthan1000bp.bed)
#  bedtools intersect -a "$infile" -b "$TEFILE" -f 0.50 -u > "$ANNOUTDIR/${base}.TEann.txt"
#done

#for infile in "$ANNOUTDIR"/*.TEann.txt; do
#  base=$(basename "$infile" .TEann.txt)
#  awk '{print $4}' "$infile" | sort | uniq -c | awk '{print $1 "\t" $2}' > "$BASEDIR/${base}_TEcounts2.bed"
#done