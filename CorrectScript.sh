#!/bin/bash
#SBATCH --job-name=BinningGenes	                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1	                                # Single task job
#SBATCH --cpus-per-task=8		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=80gb			                                # Total memory for job
#SBATCH --time=72:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/dr27977/log.%j.out
#SBATCH --error=/scratch/dr27977/log.%j.err		    # Location of standard output and error log files (replace cbergman with your myid)
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
#curl -s ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz | gunzip -c > $BASEDIR/refann1.gtf
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

#5.7.25- 5.8.25 Each row reports the number of peaks that had a certain percentage overlap with transposable elements
#module load BEDTools/2.31.0-GCC-12.3.0
#module load Homer/5.1-foss-2023a-R-4.3.2 

#PEAKS=$OUTDIR/peaks
#TE_ANN=/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/TEann_35_0.1filt.bed 
#REF_ANN=$OUTDIR/refann.gtf

#BIN_DIR=$OUTDIR/TE_overlap_bins
#mkdir -p $BIN_DIR
#SUMMARY_FILE=$BIN_DIR/bin_summary.tsv
#echo -e "sample\tbin\tcount" > $SUMMARY_FILE

#for infile in $PEAKS/*final.bed
#do
#  base=$(basename $infile _final.bed)
#  echo "Processing $base..."

#  remaining_peaks=$BIN_DIR/${base}_remaining.bed
#  cp $infile $remaining_peaks

  # Define bin ranges and labels
#  declare -a BINS=(
#    "000 0 0"
#    "010 0.0001 0.10"
#    "025 0.10 0.25"
#    "050 0.25 0.50"
#    "075 0.50 0.75"
#    "100 0.75 1.00"
#  )

#  for entry in "${BINS[@]}"
#  do
#    read -r label min max <<< "$entry"
#    bindir=$BIN_DIR/${base}_bin_${label}
#    mkdir -p $bindir
#    outfile=$bindir/${base}_TEbin_${label}.bed

#    if [[ "$label" == "000" ]]; then
#      bedtools intersect -a $remaining_peaks -b $TE_ANN -v > $outfile
#    elif [[ "$label" == "100" ]]; then
#      cp $remaining_peaks $outfile
#    else
#      bedtools intersect -a $remaining_peaks -b $TE_ANN -f $min -u > $BIN_DIR/tmp_min.bed
#      bedtools intersect -a $remaining_peaks -b $TE_ANN -f $max -v > $BIN_DIR/tmp_max.bed
#      bedtools intersect -a $BIN_DIR/tmp_min.bed -b $BIN_DIR/tmp_max.bed > $outfile
#      rm $BIN_DIR/tmp_min.bed $BIN_DIR/tmp_max.bed
#    fi

    # Remove binned peaks from remaining_peaks
#    bedtools intersect -a $remaining_peaks -b $outfile -v > $BIN_DIR/tmp_next.bed
#    mv $BIN_DIR/tmp_next.bed $remaining_peaks

    # Log number of peaks in each bin
#    count=$(wc -l < $outfile)
#    echo -e "$base\t$label\t$count" >> $SUMMARY_FILE

    # Annotate and extract genes within ±5kb of TSS
#    annotatePeaks.pl $outfile danRer11 -gtf $GTF > $bindir/${base}_TEbin_${label}.ann.txt
#    awk -F'\t' 'sqrt($10*$10) <=5000' $bindir/${base}_TEbin_${label}.ann.txt > $bindir/${base}_TEbin_${label}.within5kb.txt
#    cut -f2 $bindir/${base}_TEbin_${label}.within5kb.txt | tail -n +2 | sort | uniq > $bindir/${base}_TEbin_${label}_genes.txt
#  done

#done



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

##5.12.25 Trying again to bin genes 
#module load BEDTools
#module load Homer

#PEAKS_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks"
#TE_ANNOT="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/TEann_35_0.1filt.bed"
#REF_GTF="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/refann.gtf"
#TEBIN_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/TE_overlap_bins_by_gene"

#mkdir -p "$TEBIN_DIR"

#for peakfile in "$PEAKS_DIR"/*final.bed; do
    #base=$(basename "$peakfile" _final.bed)

    #echo "Processing $base"

    # Bin 1: 0% overlap
    #mkdir -p "$TEBIN_DIR/${base}_bin_000"
    #bedtools intersect -a "$peakfile" -b "$TE_ANNOT" -v > "$TEBIN_DIR/${base}_bin_000/${base}_TEbin_000.bed"

    # Bin 2: >0% to ≤10%
    #mkdir -p "$TEBIN_DIR/${base}_bin_010"
    #bedtools intersect -a "$peakfile" -b "$TE_ANNOT" -f 0.0001 -u > "$TEBIN_DIR/tmp_gt0.bed"
    #bedtools intersect -a "$peakfile" -b "$TE_ANNOT" -f 0.10 -u > "$TEBIN_DIR/tmp_gte10.bed"
    #bedtools intersect -a "$TEBIN_DIR/tmp_gt0.bed" -b "$TEBIN_DIR/tmp_gte10.bed" -v > "$TEBIN_DIR/${base}_bin_010/${base}_TEbin_010.bed"
    #rm "$TEBIN_DIR/tmp_gt0.bed" "$TEBIN_DIR/tmp_gte10.bed"

    # Bin 3: >10% to ≤25%
    #mkdir -p "$TEBIN_DIR/${base}_bin_025"
    #bedtools intersect -a "$peakfile" -b "$TE_ANNOT" -f 0.10 -u > "$TEBIN_DIR/tmp_gte10.bed"
    #bedtools intersect -a "$peakfile" -b "$TE_ANNOT" -f 0.25 -u > "$TEBIN_DIR/tmp_gte25.bed"
    #bedtools intersect -a "$TEBIN_DIR/tmp_gte10.bed" -b "$TEBIN_DIR/tmp_gte25.bed" -v > "$TEBIN_DIR/${base}_bin_025/${base}_TEbin_025.bed"
    #rm "$TEBIN_DIR/tmp_gte10.bed" "$TEBIN_DIR/tmp_gte25.bed"

    # Bin 4: >25% to ≤50%
    #mkdir -p "$TEBIN_DIR/${base}_bin_050"
    #bedtools intersect -a "$peakfile" -b "$TE_ANNOT" -f 0.25 -u > "$TEBIN_DIR/tmp_gte25.bed"
    #bedtools intersect -a "$peakfile" -b "$TE_ANNOT" -f 0.50 -u > "$TEBIN_DIR/tmp_gte50.bed"
    #bedtools intersect -a "$TEBIN_DIR/tmp_gte25.bed" -b "$TEBIN_DIR/tmp_gte50.bed" -v > "$TEBIN_DIR/${base}_bin_050/${base}_TEbin_050.bed"
    #rm "$TEBIN_DIR/tmp_gte25.bed" "$TEBIN_DIR/tmp_gte50.bed"

    ## Bin 5: >50% to ≤75%
    #mkdir -p "$TEBIN_DIR/${base}_bin_075"
    #bedtools intersect -a "$peakfile" -b "$TE_ANNOT" -f 0.50 -u > "$TEBIN_DIR/tmp_gte50.bed"
    #bedtools intersect -a "$peakfile" -b "$TE_ANNOT" -f 0.75 -u > "$TEBIN_DIR/tmp_gte75.bed"
    #bedtools intersect -a "$TEBIN_DIR/tmp_gte50.bed" -b "$TEBIN_DIR/tmp_gte75.bed" -v > "$TEBIN_DIR/${base}_bin_075/${base}_TEbin_075.bed"
    #rm "$TEBIN_DIR/tmp_gte50.bed" "$TEBIN_DIR/tmp_gte75.bed"

    # Bin 6: >75% to ≤100%
    #mkdir -p "$TEBIN_DIR/${base}_bin_100"
    #bedtools intersect -a "$peakfile" -b "$TE_ANNOT" -f 0.75 -u > "$TEBIN_DIR/${base}_bin_100/${base}_TEbin_100.bed"

    ## Annotate each bin and extract gene names within ±5kb of TSS
#    for binfile in "$TEBIN_DIR/${base}_bin_"*/${base}_TEbin_*.bed; do
#        bindir=$(dirname "$binfile")
#        binbase=$(basename "$binfile" .bed)

#        annotatePeaks.pl "$binfile" danRer11 -gtf "$REF_GTF" > "$bindir/${binbase}.ann.txt"
#        awk -F'\t' 'NR==1 {print} NR>1 && sqrt($10*$10) <= 5000' "$bindir/${binbase}.ann.txt" > "$bindir/${binbase}.within5kb.txt"

#        cut -f2 "$bindir/${binbase}.within5kb.txt" | tail -n +2 | sort | uniq > "$bindir/${binbase}_genes.txt"
#    done

#done

##5.12.25
#module load BEDTools
#module load Homer

# === CONFIGURATION ===
#PEAKS_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks"
#TE_ANNOT="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/TEann_35_0.1filt.bed"
#REF_GTF="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/refann.gtf"
#TEBIN_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/TE_overlap_bins_by_gene2"
#mkdir -p "$TEBIN_DIR"

# === BIN THRESHOLDS ===
#declare -a BINS=("000" "010" "025" "050" "075" "100")
#declare -a F_LOWER=(0.00 0.0001 0.10 0.25 0.50 0.75)
#declare -a F_UPPER=(0.00 0.10   0.25 0.50 0.75 1.01)

#for peakfile in "$PEAKS_DIR"/*final.bed; do
#    base=$(basename "$peakfile" _final.bed)
#    echo "Processing $base..."

#    remaining="$peakfile"

#    for i in "${!BINS[@]}"; do
#        bin="${BINS[$i]}"
#        flo="${F_LOWER[$i]}"
#        fhi="${F_UPPER[$i]}"
#        bindir="$TEBIN_DIR/${base}_bin_${bin}"
#        mkdir -p "$bindir"
#        outfile="$bindir/${base}_TEbin_${bin}.bed"

        # Get peaks within this overlap fraction
#        bedtools intersect -a "$remaining" -b "$TE_ANNOT" -f "$flo" -F 0.01 -u > "$outfile.tmp"

#        if [[ "$bin" == "000" ]]; then
#            bedtools intersect -a "$remaining" -b "$TE_ANNOT" -v > "$outfile"
#        else
#            bedtools intersect -a "$outfile.tmp" -b "$TE_ANNOT" -f "$fhi" -F 0.01 -v > "$outfile"
#        fi

        # Subtract this bin’s peaks from the remaining set
#        remaining_tmp="$bindir/${base}_remaining.bed"
#        bedtools subtract -a "$remaining" -b "$outfile" > "$remaining_tmp"
#        mv "$remaining_tmp" "$remaining"
#        rm -f "$outfile.tmp"
#    done

    # === ANNOTATE AND EXTRACT VALID GENE SYMBOLS WITHIN ±5kb ===
#    for binfile in "$TEBIN_DIR/${base}_bin_"*/${base}_TEbin_*.bed; do
#        bindir=$(dirname "$binfile")
#        binbase=$(basename "$binfile" .bed)

#        annotatePeaks.pl "$binfile" danRer11 -gtf "$REF_GTF" > "$bindir/${binbase}.ann.txt"
#        awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="Gene Name") col=i} NR>1 && col && sqrt($10*$10) <= 5000 && $col != "NA" && $col != "." {print $col}' "$bindir/${binbase}.ann.txt" | sort | uniq > "$bindir/${binbase}_genes.txt"
#    done
#done

#TEBIN_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/TE_overlap_bins_by_gene2"
#echo -e "Sample\tBin\tGene_Count"

# Loop through all _genes.txt files
#find "$TEBIN_DIR" -name "*_genes.txt" | sort | while read gene_file; do
#    count=$(wc -l < "$gene_file")
#    base=$(basename "$gene_file" _genes.txt)
#    folder=$(basename "$(dirname "$gene_file")")
#    sample=${folder%%_bin_*}
#    bin=${folder##*_bin_}
#    echo -e "${sample}\t${bin}\t${count}"
#done

##5.13.25
#module load BEDTools
#module load Homer

#PEAKS_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks"
#TE_BED="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/TEann_35_0.1filt.bed"
#GTF="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/refann.gtf"
#GENIC_TE_BIN_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/Genic_TE_overlap_bins"
#LOG_DIR="/scratch/dr27977/logs"

#mkdir -p "$GENIC_TE_BIN_DIR"
#mkdir -p "$LOG_DIR"

#for peakfile in "$PEAKS_DIR"/*final.bed; do
#    base=$(basename "$peakfile" _final.bed)
#    echo "Processing $base"

#    # Annotate peaks and filter to genic peaks within 5kb of TSS
#    annotatePeaks.pl "$peakfile" danRer11 -gtf "$GTF" > "$GENIC_TE_BIN_DIR/${base}.ann.txt"
#    awk -F'\t' 'NR==1 || sqrt($10*$10) <= 5000' "$GENIC_TE_BIN_DIR/${base}.ann.txt" > "$GENIC_TE_BIN_DIR/${base}_within5kb.txt"
#    awk 'NR > 1 {print $2"\t"$3"\t"$4"\t"$1}' "$GENIC_TE_BIN_DIR/${base}_within5kb.txt" > "$GENIC_TE_BIN_DIR/${base}_genic.bed"

    # Bin 1: 0% TE overlap
#    mkdir -p "$GENIC_TE_BIN_DIR/${base}_bin_000"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/${base}_genic.bed" -b "$TE_BED" -v > "$GENIC_TE_BIN_DIR/${base}_bin_000/${base}_TEbin_000.bed"

    # Bin 2: >0% to ≤10% TE overlap
#    mkdir -p "$GENIC_TE_BIN_DIR/${base}_bin_010"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/${base}_genic.bed" -b "$TE_BED" -f 0.0001 -u > "$GENIC_TE_BIN_DIR/tmp_gt0.bed"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/${base}_genic.bed" -b "$TE_BED" -f 0.10 -u > "$GENIC_TE_BIN_DIR/tmp_gte10.bed"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/tmp_gt0.bed" -b "$GENIC_TE_BIN_DIR/tmp_gte10.bed" -v > "$GENIC_TE_BIN_DIR/${base}_bin_010/${base}_TEbin_010.bed"
#    rm "$GENIC_TE_BIN_DIR/tmp_gt0.bed" "$GENIC_TE_BIN_DIR/tmp_gte10.bed"

    # Bin 3: >10% to ≤25% TE overlap
#    mkdir -p "$GENIC_TE_BIN_DIR/${base}_bin_025"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/${base}_genic.bed" -b "$TE_BED" -f 0.10 -u > "$GENIC_TE_BIN_DIR/tmp_gte10.bed"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/${base}_genic.bed" -b "$TE_BED" -f 0.25 -u > "$GENIC_TE_BIN_DIR/tmp_gte25.bed"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/tmp_gte10.bed" -b "$GENIC_TE_BIN_DIR/tmp_gte25.bed" -v > "$GENIC_TE_BIN_DIR/${base}_bin_025/${base}_TEbin_025.bed"
#    rm "$GENIC_TE_BIN_DIR/tmp_gte10.bed" "$GENIC_TE_BIN_DIR/tmp_gte25.bed"

    # Bin 4: >25% to ≤50% TE overlap
#    mkdir -p "$GENIC_TE_BIN_DIR/${base}_bin_050"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/${base}_genic.bed" -b "$TE_BED" -f 0.25 -u > "$GENIC_TE_BIN_DIR/tmp_gte25.bed"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/${base}_genic.bed" -b "$TE_BED" -f 0.50 -u > "$GENIC_TE_BIN_DIR/tmp_gte50.bed"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/tmp_gte25.bed" -b "$GENIC_TE_BIN_DIR/tmp_gte50.bed" -v > "$GENIC_TE_BIN_DIR/${base}_bin_050/${base}_TEbin_050.bed"
#    rm "$GENIC_TE_BIN_DIR/tmp_gte25.bed" "$GENIC_TE_BIN_DIR/tmp_gte50.bed"

    # Bin 5: >50% to ≤75% TE overlap
#    mkdir -p "$GENIC_TE_BIN_DIR/${base}_bin_075"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/${base}_genic.bed" -b "$TE_BED" -f 0.50 -u > "$GENIC_TE_BIN_DIR/tmp_gte50.bed"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/${base}_genic.bed" -b "$TE_BED" -f 0.75 -u > "$GENIC_TE_BIN_DIR/tmp_gte75.bed"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/tmp_gte50.bed" -b "$GENIC_TE_BIN_DIR/tmp_gte75.bed" -v > "$GENIC_TE_BIN_DIR/${base}_bin_075/${base}_TEbin_075.bed"
#    rm "$GENIC_TE_BIN_DIR/tmp_gte50.bed" "$GENIC_TE_BIN_DIR/tmp_gte75.bed"

    # Bin 6: >75% to ≤100% TE overlap
#    mkdir -p "$GENIC_TE_BIN_DIR/${base}_bin_100"
#    bedtools intersect -a "$GENIC_TE_BIN_DIR/${base}_genic.bed" -b "$TE_BED" -f 0.75 -u > "$GENIC_TE_BIN_DIR/${base}_bin_100/${base}_TEbin_100.bed"

    # Annotate and extract gene names
#    for bedfile in "$GENIC_TE_BIN_DIR/${base}_bin_"*/${base}_TEbin_*.bed; do
#        dir=$(dirname "$bedfile")
#        name=$(basename "$bedfile" .bed)

#        annotatePeaks.pl "$bedfile" danRer11 -gtf "$GTF" > "$dir/${name}.ann.txt"
#        awk -F'\t' 'NR==1 || sqrt($10*$10) <= 5000' "$dir/${name}.ann.txt" > "$dir/${name}.within5kb.txt"
#        cut -f2 "$dir/${name}.within5kb.txt" | tail -n +2 | sort | uniq > "$dir/${name}_genes.txt"
#    done

#    echo "Finished $base"
#done

#echo "Binning complete. Results saved to $GENIC_TE_BIN_DIR"

#module load BEDTools
#module load Homer

#PEAKS_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks"
#TE_ANNOT="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/TEann_35_0.1filt.bed"
#REF_GTF="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/refann.gtf"
#OUT_BASE="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/H3K9me3_TSS_TE_categories_with_genes"

#mkdir -p "$OUT_BASE"

#for peakfile in "$PEAKS_DIR"/*final.bed; do
#    base=$(basename "$peakfile" _final.bed)
#    echo "Processing $base"

    # Annotate and filter to within ±5kb of TSS
#    annotatePeaks.pl "$peakfile" danRer11 -gtf "$REF_GTF" > "$OUT_BASE/${base}.ann.txt"
#    awk -F'\t' 'NR==1 || ($8 ~ /exon|intron/ && sqrt($10*$10) <= 5000)' "$OUT_BASE/${base}.ann.txt" > "$OUT_BASE/${base}_within5kb.txt"

    # Split into categories
#    awk -F'\t' '
#        NR==1 {next}
#        $8 ~ /exon/ && $8 !~ /intron/ {print $2 "\t" $3 "\t" $4 > "'"$OUT_BASE/${base}_exon_only.bed"'";}
#        $8 ~ /intron/ && $8 !~ /exon/ {print $2 "\t" $3 "\t" $4 > "'"$OUT_BASE/${base}_intron_only.bed"'";}
#        $8 ~ /intron/ && $8 ~ /exon/  {print $2 "\t" $3 "\t" $4 > "'"$OUT_BASE/${base}_both.bed"'";}
#    ' "$OUT_BASE/${base}_within5kb.txt"

    # Process each category
#    for cat in exon_only intron_only both; do
#        mkdir -p "$OUT_BASE/${base}_${cat}"

#        sort -k1,1 -k2,2n "$OUT_BASE/${base}_${cat}.bed" > "$OUT_BASE/${base}_${cat}/input.bed"

        # Bin 000 (0% TE overlap)
#        bedtools intersect -a "$OUT_BASE/${base}_${cat}/input.bed" -b "$TE_ANNOT" -v \
#          > "$OUT_BASE/${base}_${cat}/${cat}_bin_000.bed"

        # Bins: ≤10, ≤25, ≤50, ≤75
#        for threshold in 0.10 0.25 0.50 0.75; do
#            perc=$(echo "$threshold" | awk '{printf "%03d", $1*100}')
#            binname="${cat}_bin_${perc}"
#            bedtools intersect -a "$OUT_BASE/${base}_${cat}/input.bed" -b "$TE_ANNOT" -f "$threshold" -u \
#              > "$OUT_BASE/${base}_${cat}/${binname}.bed"
#        done

        # Bin 100 = everything else not in previous bins
#        cat "$OUT_BASE/${base}_${cat}/${cat}_bin_"*.bed 2>/dev/null | sort -k1,1 -k2,2n | uniq \
#          > "$OUT_BASE/${base}_${cat}/tmp_all_bins_covered.bed"
#        bedtools subtract -a "$OUT_BASE/${base}_${cat}/input.bed" -b "$OUT_BASE/${base}_${cat}/tmp_all_bins_covered.bed" \
#          > "$OUT_BASE/${base}_${cat}/${cat}_bin_100.bed"
#        rm "$OUT_BASE/${base}_${cat}/tmp_all_bins_covered.bed"

        # Annotate each bin and extract gene symbols
#        for binfile in "$OUT_BASE/${base}_${cat}/${cat}_bin_"*.bed; do
#            binbase=$(basename "$binfile" .bed)
#            annotatePeaks.pl "$binfile" danRer11 -gtf "$REF_GTF" > "$OUT_BASE/${base}_${cat}/${binbase}.ann.txt"
#            awk -F'\t' 'NR > 1 && $2 != "NA" {print $2}' "$OUT_BASE/${base}_${cat}/${binbase}.ann.txt" | sort | uniq \
#              > "$OUT_BASE/${base}_${cat}/${binbase}_genes.txt"
#        done
#    done
#done

##5.14.25 to make it extract gene names 
#module load mygene
#module load SciPy-bundle

#IN_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/H3K9me3_TSS_TE_categories_with_genes"
#OUT_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/gene_symbol_lookup_output"

#mkdir -p "$OUT_DIR"

#python <<EOF
#import os
#import pandas as pd
#from mygene import MyGeneInfo

#mg = MyGeneInfo()
#in_dir = "$IN_DIR"
#out_dir = "$OUT_DIR"

#for file in os.listdir(in_dir):
#    if file.endswith("_within5kb.txt"):
#        in_path = os.path.join(in_dir, file)
#        out_path = os.path.join(out_dir, file.replace(".txt", "_symbols.txt"))
#        print(f"Processing {file}")

        # Extract ENSDART/ENSDARG IDs from columns 11 and 12
#        ens_ids = set()
#        with open(in_path) as f:
#            header = f.readline()
#            for line in f:
#                fields = line.strip().split("\t")
#                if len(fields) > 12:
#                    for i in [10, 11]:  # promoter and gene IDs
#                        if fields[i].startswith("ENSDAR"):
#                            ens_ids.add(fields[i])

#        if not ens_ids:
#            print(f"No valid Ensembl IDs in {file}")
#            continue

#        results = mg.querymany(list(ens_ids), scopes="ensembl.gene", fields="symbol", species="zebrafish", as_dataframe=True)

#        if not results.empty and "symbol" in results.columns:
#            df = results[["symbol"]].reset_index().rename(columns={"query": "Ensembl_ID", "symbol": "Gene_Symbol"})
#            df.to_csv(out_path, sep="\t", index=False)
#        else:
#            print(f"No mapping found for {file}")
#EOF

## make a summary table

##5.27.25 Trying again to make summary tables
#module load BEDTools
#module load Homer

#BASEDIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published"
#PEAKS_DIR="$BASEDIR/peaks"
#TE_BED="$PEAKS_DIR/TEann_35_0.1filt.sorted.bed"
#GTF="$BASEDIR/refann.gtf"
#OUTDIR="$BASEDIR/H3K9me3_summary_tables"
#TMPDIR="$OUTDIR/tmp"
#GENOME="danRer11"
#CHRLEN="$BASEDIR/genome/chrNameLength.txt"

#[ -f "$GTF" ] || { echo "ERROR: GTF not found: $GTF"; exit 1; }
#[ -f "$CHRLEN" ] || { echo "ERROR: chrNameLength.txt not found: $CHRLEN"; exit 1; }

#mkdir -p "$OUTDIR" "$TMPDIR"

#if [ ! -f "$TE_BED" ]; then
#    echo "Sorting TE BED file..."
#    bedtools sort -i "$PEAKS_DIR/TEann_35_0.1filt.bed" -g "$CHRLEN" > "$TE_BED"
#fi

#shopt -s nullglob
#peakfiles=("$PEAKS_DIR"/*final.bed)

#if [ ${#peakfiles[@]} -eq 0 ]; then
#    echo "ERROR: No peak files found."
#    exit 1
#fi

#for peakfile in "${peakfiles[@]}"; do
#    base=$(basename "$peakfile" _final.bed)
#    echo "Processing $base..."

#    for window in 1000 5000; do
#        echo "  TSS window: ${window}bp"

#        annfile="$TMPDIR/${base}_ann_${window}.txt"
#        filtered="$TMPDIR/${base}_filtered_${window}.bed"
#        covfile="$TMPDIR/${base}_cov_${window}.txt"
#        outtable="$OUTDIR/${base}_TSS${window}bp_TE_table.tsv"

#        echo "    Annotating peaks with HOMER..."
#        annotatePeaks.pl "$peakfile" $GENOME -gtf "$GTF" > "$annfile"
#        if [ ! -s "$annfile" ]; then
#            echo "    ERROR: Annotation failed: $annfile"
#            continue
#        fi

#        echo "    Filtering and formatting..."
#        awk -v W=$window 'BEGIN{OFS="\t"}
#        NR > 1 && sqrt($10*$10) <= W {
#            chr=$2; start=$3; end=$4; peakID=$1;
#            category="unannotated";
#            if ($8 ~ /exon/ && $8 !~ /intron/) category="exon_only";
#            else if ($8 ~ /intron/ && $8 !~ /exon/) category="intron_only";
#            else if ($8 ~ /exon/ && $8 ~ /intron/) category="both";
#            else if ($8 ~ /firstExon/) category="first_exon";

#            print chr, start, end, peakID, category, $10, $11, $12
#        }' "$annfile" > "$filtered.unsorted"

#        bedtools sort -i "$filtered.unsorted" -g "$CHRLEN" > "$filtered"
#        rm "$filtered.unsorted"

#        echo "    Calculating TE overlap..."
#        bedtools coverage -a "$filtered" -b "$TE_BED" -sorted -f 0.0001 > "$covfile"
#        if [ ! -s "$covfile" ]; then
#            echo "    ERROR: TE coverage failed: $covfile"
#            continue
#        fi

#        echo "    Writing output table..."
#        awk -v time="$base" -v win="$window" '
#        BEGIN { OFS="\t"; print "gene_id","gene_name","peak_id","TSS_dist","category","TE_overlap_pct","TE_bin","timepoint" }
#        {
#            gene_id = $7;
#            gene_name = $8;
#            peak_id = $4;
#            tss_dist = $6;
#            category = $5;
#            overlap_pct = $NF;
#            bin = "100%";

#            if (overlap_pct == 0) bin = "0%";
#            else if (overlap_pct <= 10) bin = "<=10%";
#            else if (overlap_pct <= 25) bin = "<=25%";
#            else if (overlap_pct <= 50) bin = "<=50%";
#            else if (overlap_pct <= 75) bin = "<=75%";

#            print gene_id, gene_name, peak_id, tss_dist, category, overlap_pct, bin, time
#        }' "$covfile" > "$outtable"

#        echo "    Output written: $outtable"
#    done
#done

#echo "All timepoints processed."

##5.27.25 Get gene names
#GTF="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/refann.gtf"
#INPUT_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/H3K9me3_summary_tables"
#OUTPUT_DIR="$INPUT_DIR/with_symbols"
#LOOKUP="$OUTPUT_DIR/zebrafish_ensid_to_symbol.tsv"
#mkdir -p "$OUTPUT_DIR"

# === STEP 1: Extract gene symbol mappings from GTF ===
#echo "Extracting gene_id and transcript_id to gene_name mappings from GTF..."

#awk -F'\t' '
#$3 == "gene" || $3 == "transcript" {
#    match($9, /gene_id "([^"]+)"/, gid);
#    match($9, /transcript_id "([^"]+)"/, tid);
#    match($9, /gene_name "([^"]+)"/, gname);
#    if (gid[1] && gname[1]) print gid[1] "\t" gname[1];
#    if (tid[1] && gname[1]) print tid[1] "\t" gname[1];
#}' "$GTF" | sort -u > "$LOOKUP"

#echo "Lookup table saved to: $LOOKUP"

# === STEP 2: Annotate all TE_table files with gene symbols ===
#echo "Annotating TE tables with gene symbols..."

#for file in "$INPUT_DIR"/*_TE_table.tsv; do
#    [ -f "$file" ] || continue
#    base=$(basename "$file" .tsv)
#    output="$OUTPUT_DIR/${base}_with_symbols.tsv"

#    awk -F'\t' -v OFS='\t' -v mapfile="$LOOKUP" '
#    BEGIN {
#        while ((getline < mapfile) > 0) {
#            id_to_name[$1] = $2;
#        }
#    }
#    NR == 1 {
#        print $0, "gene_symbol";
#        next;
#    }
#    {
#        symbol = ($1 in id_to_name) ? id_to_name[$1] : "NA";
#        print $0, symbol;
#    }' "$file" > "$output"

#    echo "Annotated: $output"
#done

#echo "All TE tables processed with gene symbols added."
#BASE_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/H3K9me3_summary_tables/with_symbols"
#OUTPUT_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/H3K9me3_summary_tables"
#OUTPUT_1KB="${OUTPUT_DIR}/gene_peak_summaries_1kb.txt"
#OUTPUT_5KB="${OUTPUT_DIR}/gene_peak_summaries_5kb.txt"

#rm -f "$OUTPUT_1KB" "$OUTPUT_5KB"

#for file in "${BASE_DIR}"/*_TSS1000bp_TE_table_with_symbols.tsv "${BASE_DIR}"/*_TSS5000bp_TE_table_with_symbols.tsv; do
#    [ -f "$file" ] || continue

#    filename=$(basename "$file")
#    if [[ "$filename" == *_TSS1000bp_* ]]; then
#        OUTPUT="$OUTPUT_1KB"
#    elif [[ "$filename" == *_TSS5000bp_* ]]; then
#        OUTPUT="$OUTPUT_5KB"
#    else
#        continue
#    fi

#    awk -F'\t' -v OFS='\t' -v out="$OUTPUT" '
#    BEGIN {
#        categories["intron_only"]
#        categories["exon_only"]
#        categories["both"]
#        categories["first_exon"]
#        categories["unannotated"]
#    }

#    NR == 1 {
#        for (i = 1; i <= NF; i++) header[$i] = i
#        next
#    }

#    {
#        raw_gid = $header["gene_id"]
#        gsym = $header["gene_symbol"]
#        tp = $header["timepoint"]
#        cat = $header["category"]
#        tss = $header["TSS_dist"]
#        te = $header["TE_overlap_pct"]

#        gsub(/[()]/, "", raw_gid)
#        split(raw_gid, a, " ")
#        gid = a[1]

#        if (gid == "" || tp == "" || cat == "") next

#        key = gid "|" gsym "|" tp

#        count[key]++
#        cat_count[key, cat]++
#        if (te != "" && te ~ /^[0-9.]+$/) {
#            te_min[key] = (key in te_min) ? (te < te_min[key] ? te : te_min[key]) : te
#            te_max[key] = (key in te_max) ? (te > te_max[key] ? te : te_max[key]) : te
#        }

#        if (tss ~ /^-?[0-9.]+$/ && tss + 0 >= -1000 && tss + 0 <= 1000) {
#            tss_peak[key]++
#            tss_pos[key] = (tss_pos[key] ? tss_pos[key] ", " : "") tss
#        }
#    }

#    END {
#        for (k in count) {
#            split(k, a, "|")
#            gid = a[1]; gsym = a[2]; tp = a[3]

#            print "Transcript ID: " gid >> out
#            print "Gene Symbol: " (gsym == "" ? "NA" : gsym) >> out
#            print "Timepoint: " tp >> out

#            if (tss_peak[k] > 0) {
#                print "TSS peak: Yes (" tss_peak[k] " peak[s] at: " tss_pos[k] ")" >> out
#            } else {
#                print "TSS peak: No" >> out
#            }

#            print "H3K9me3 peaks:" >> out
#            for (c in categories) {
#                val = (cat_count[k, c] ? cat_count[k, c] : 0)
#                print "  " val " " c " peaks" >> out
#            }

#            if (k in te_min && k in te_max) {
#                printf("TE overlap range: %.2f - %.2f\n", te_min[k], te_max[k]) >> out
#            } else {
#                print "TE overlap range: NA" >> out
#            }

#            print "------------------------------------------------------------" >> out
#        }
#    }' "$file"
#done

#echo "Finished summarizing 1kb and 5kb gene peak tables."

##5.28.25 Trying to get bins of genes in excel files/ 5.29.25
#module load BEDTools/2.31.0-GCC-12.3.0
#module load Homer/5.1-foss-2023a-R-4.3.2

#BED_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks"
#GTF="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/refann1.gtf"
#OUT_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/finalsummaries5"
#TMP="${OUT_DIR}/tmp"

#mkdir -p "$TMP"
#mkdir -p "$OUT_DIR"

# Loop through each BED file
#for bedfile in "$BED_DIR"/*_K9_final.bed; do
#    base=$(basename "$bedfile" _K9_final.bed)

#    for WIN in 1000 5000; do
#        window_label="$( [ "$WIN" -eq 1000 ] && echo "1kb" || echo "5kb")"
#        trimmed_bed="$TMP/${base}_${window_label}_trimmed.bed"
#        annot_file="$TMP/${base}_${window_label}_annot.txt"
#        output_tsv="$OUT_DIR/${base}_K9_TSS${WIN}bp_TE_table_with_symbols_summary.tsv"

        # Trim to first 6 columns
#        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' "$bedfile" > "$trimmed_bed"

        # Annotate with HOMER
#        annotatePeaks.pl "$trimmed_bed" danRer11 -gtf "$GTF" -size $WIN > "$annot_file"

        # Run external Python script
#        python3 compute_TE_overlap.py "$annot_file" "$trimmed_bed" "$output_tsv" "$base" "$window_label"
#    done
#done

#echo "All done. Output saved to: $OUT_DIR"

##6.3.25 From Scratch,Doing what Katie and Ashley did and then working off of that script
#module load BEDTools/2.31.0-GCC-12.3.0
#module load Homer/5.1-foss-2023a-R-4.3.2 

#BASE_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published"
#MASKANN_DIR="$BASE_DIR/peaks/ann"
#TE_BED="$BASE_DIR/peaks/TEann_35_0.1filt.bed"
#GTF="$BASE_DIR/refann.gtf"

# Create output directories
#OUT_1KB="$BASE_DIR/peaks/filtered_peaks_1kb_noTE"
#OUT_5KB="$BASE_DIR/peaks/filtered_peaks_5kb_noTE"
#mkdir -p "$OUT_1KB" "$OUT_5KB"

############################################
# Step 1: Filter for distance to TSS (1kb and 5kb)
############################################

# 1kb filtering
#for file in $MASKANN_DIR/*maskann.txt
#do
#  base=$(basename "$file" _maskann.txt)
#  awk -F'\t' 'NR==1 || sqrt($10*$10) <= 1000' "$file" > "$OUT_1KB/${base}.1000bp_ann.txt"
#done

# 5kb filtering
#for file in $MASKANN_DIR/*maskann.txt
#do
#  base=$(basename "$file" _maskann.txt)
#  awk -F'\t' 'NR==1 || sqrt($10*$10) <= 5000' "$file" > "$OUT_5KB/${base}.5000bp_ann.txt"
#done

############################################
# Step 2: Convert filtered .txt to BED format
############################################

# 1kb BED
#for file in $OUT_1KB/*.1000bp_ann.txt
#do
#  base=$(basename "$file" .1000bp_ann.txt)
#  awk -F'\t' 'NR>1 {OFS="\t"; print $2, $3, $4, $1}' "$file" > "$OUT_1KB/${base}.1000bp.bed"
#done

# 5kb BED
#for file in $OUT_5KB/*.5000bp_ann.txt
#do
#  base=$(basename "$file" .5000bp_ann.txt)
#  awk -F'\t' 'NR>1 {OFS="\t"; print $2, $3, $4, $1}' "$file" > "$OUT_5KB/${base}.5000bp.bed"
#done

############################################
# Step 3: Remove peaks overlapping TEs by ≥50%
############################################

# 1kb noTE
#for bed in $OUT_1KB/*.1000bp.bed
#do
#  base=$(basename "$bed" .1000bp.bed)
#  bedtools intersect -a "$bed" -b "$TE_BED" -f 0.50 -v > "$OUT_1KB/${base}.1000bp_noTE.bed"
#done

# 5kb noTE
#for bed in $OUT_5KB/*.5000bp.bed
#do
#  base=$(basename "$bed" .5000bp.bed)
#  bedtools intersect -a "$bed" -b "$TE_BED" -f 0.50 -v > "$OUT_5KB/${base}.5000bp_noTE.bed"
#done

############################################
# Step 4: (Optional) Re-annotate filtered noTE peaks with HOMER
############################################

#for bed in $OUT_1KB/*.1000bp_noTE.bed
#do
#  base=$(basename "$bed" .1000bp_noTE.bed)
#  annotatePeaks.pl "$bed" danRer11 -gtf "$GTF" > "$OUT_1KB/${base}.1000bp_noTE.ann.txt"
#done

#for bed in $OUT_5KB/*.5000bp_noTE.bed
#do
#  base=$(basename "$bed" .5000bp_noTE.bed)
#  annotatePeaks.pl "$bed" danRer11 -gtf "$GTF" > "$OUT_5KB/${base}.5000bp_noTE.ann.txt"
#done

#done

##6.4.24 Making files that show peaks within 1kb and within 5kb of a TSS (SUCCESSFUL)
#module load Homer/5.1-foss-2023a-R-4.3.2
#curl -s ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz | gunzip -c > $OUTDIR/refann.gtf
#mkdir -p $OUTDIR/peaksnew/ann


# for infile in /scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/*final.bed
# do
#   base=$( basename ${infile} _final.bed)
#   annotatePeaks.pl $infile danRer11 -gtf $OUTDIR/refann.gtf > $OUTDIR/peaksnew/ann/${base}.maskann.txt
# done

##Filter peaks within 1kb of TSS
# for infile in $OUTDIR/peaksnew/ann/*maskann.txt
# do
#   base=$(basename ${infile} .maskann.txt)
#   awk -F'\t' 'sqrt($10*$10) <=1000' $infile > $OUTDIR/peaksnew/ann/${base}.1000bp_ann.txt
# done
## Filter peaks greater than 1Kb of TSS
# for infile in $OUTDIR/peaksnew/ann/*maskann.txt
# do
#   base=$(basename ${infile} .maskann.txt)
#   awk -F'\t' 'sqrt($10*$10) >=1000' $infile | awk '{print $2 "\t" $3 "\t" $4 }' > $OUTDIR/peaksnew/ann/${base}.MOREthan1000bp.bed
# done

 ##Filter peaks within 5kb of TSS
# for infile in $OUTDIR/peaksnew/ann/*maskann.txt
# do
#   base=$(basename ${infile} .maskann.txt)
#   awk -F'\t' 'sqrt($10*$10) <=5000' $infile > $OUTDIR/peaksnew/ann/${base}.5000bp_ann.txt
# done
## Filter peaks greater than 5Kb of TSS
# for infile in $OUTDIR/peaksnew/ann/*maskann.txt
# do
#   base=$(basename ${infile} .maskann.txt)
#   awk -F'\t' 'sqrt($10*$10) >=5000' $infile | awk '{print $2 "\t" $3 "\t" $4 }' > $OUTDIR/peaksnew/ann/${base}.MOREthan5000bp.bed
# done

##6.5.25 Trying to add in TE component (Successful I think)
#module load BEDTools/2.31.0-GCC-12.3.0
# ==== SET PATHS ====
#ANN_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaksnew/ann"
#TE_BED="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/TEann_35_0.1filt.bed"
#GTF="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/refann.gtf"
#OUT_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaksnew/TE/TEwithorwithoutoverlap"
#mkdir -p $OUT_DIR

# ==== STEP 0: Extract transcript ID → gene symbol from GTF ====
#MAP_TSV="$OUT_DIR/zebrafish_ensid_to_symbol.tsv"

#awk -F'\t' '$3 == "transcript" && $9 ~ /transcript_id/ && $9 ~ /gene_name/ {
#    match($9, /transcript_id "([^"]+)"/, tid)
#    match($9, /gene_name "([^"]+)"/, gname)
#    if (tid[1] && gname[1]) {
#        print tid[1] "\t" gname[1]
#    }
#}' "$GTF" | sort -u > "$MAP_TSV"

# ==== PROCESS EACH 1kb & 5kb PEAK FILE ====
#for annfile in $ANN_DIR/*.1000bp_ann.txt $ANN_DIR/*.5000bp_ann.txt; do
#    base=$(basename $annfile .txt)
#    echo "Processing $base..."

    # Step 1: Extract peakID → transcript ID from HOMER annotation
#    awk 'NR > 1 {print "peak"NR, $12}' "$annfile" > "$OUT_DIR/${base}_peak_to_transcript.tsv"

    # Step 2: Join with transcript ID → gene symbol from GTF
#    join -1 2 -2 1 <(sort -k2,2 "$OUT_DIR/${base}_peak_to_transcript.tsv") <(sort -k1,1 "$MAP_TSV") \
#        | awk '{print $2 "\t" $1 "\t" $3}' > "$OUT_DIR/${base}_peak_gene_map.tsv"
    # Format: peakID    transcriptID    geneSymbol

    # Step 3: Convert HOMER to BED
#    awk 'NR > 1 {OFS="\t"; print $2, $3, $4, "peak"NR}' "$annfile" > "$OUT_DIR/${base}.bed"

    # Step 4: Intersect with TE BED (include non-overlapping peaks)
#    bedtools intersect -a "$OUT_DIR/${base}.bed" -b "$TE_BED" -wao > "$OUT_DIR/${base}_TE_overlap_raw.txt"

    # Step 5: Annotate intersect output with GeneID and GeneSymbol
#    awk '
#    BEGIN {
#        while ((getline < "'"$OUT_DIR/${base}_peak_gene_map.tsv"'") > 0) {
#            map[$1] = $2 "\t" $3
#        }
#    }
#    {
#        key = $4
#        gene = (key in map) ? map[key] : "NA\tNA"
#        print $0 "\t" gene
#    }' "$OUT_DIR/${base}_TE_overlap_raw.txt" > "$OUT_DIR/${base}_TE_overlap_with_genes.txt"

    # Step 6: Output WITH TE overlap
#    {
#        echo -e "Chr\tStart\tEnd\tPeakID\tTE_Chr\tTE_Start\tTE_End\tTE_Name\tTE_Score\tOverlap_bp\tGeneID\tGeneSymbol"
#        awk '$10 > 0' "$OUT_DIR/${base}_TE_overlap_with_genes.txt"
#    } > "$OUT_DIR/${base}_WITH_TE_overlap.txt"

    # Step 7: Output NO TE overlap
#    {
#        echo -e "Chr\tStart\tEnd\tPeakID\tTE_Chr\tTE_Start\tTE_End\tTE_Name\tTE_Score\tOverlap_bp\tGeneID\tGeneSymbol"
#        awk '$10 == 0' "$OUT_DIR/${base}_TE_overlap_with_genes.txt"
#    } > "$OUT_DIR/${base}_NO_TE_overlap.txt"
#
#done

#echo "Done. Output saved to: $OUT_DIR"

##6.5.25 trying to get summary tables 
#module load BEDTools/2.31.0-GCC-12.3.0

# ==== SET PATHS ====
#ANN_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaksnew/ann"
#TE_BED="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/TEann_35_0.1filt.bed"
#OUT_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaksnew/SummaryTables"
#mkdir -p "$OUT_DIR"

# ==== PROCESS EACH PEAK FILE ====
#for annfile in $ANN_DIR/*.1000bp_ann.txt $ANN_DIR/*.5000bp_ann.txt; do
#    base=$(basename "$annfile" .txt)
#    echo "Processing $base"

    # Generate temp BED with PeakID
#    awk 'NR>1 {OFS="\t"; print $2, $3, $4, "peak"NR}' "$annfile" > ${base}.bed

    # Extract PeakID to GeneID and GeneSymbol, and exon/intron classification
#    awk -F'\t' '
#    NR > 1 {
#        peak = "peak"NR
#        gene_id = $12
#        gene_symbol = ($17 == "") ? "NA" : $17
#        exon = ($9 ~ /exon/) ? 1 : 0
#        intron = ($9 ~ /intron/) ? 1 : 0
#        class = (exon && !intron) ? "exon_only" :
#                (!exon && intron) ? "intron_only" :
#                (exon && intron) ? "exon+intron" : "unclassified"
#        print peak "\t" gene_id "\t" gene_symbol "\t" class
#    }' "$annfile" > ${base}_peak_gene_class.tsv

    # Intersect with TE BED (include non-overlapping)
#    bedtools intersect -a ${base}.bed -b "$TE_BED" -wao > ${base}_TEraw.txt

    # Summarize TE overlap: percent, bin, and TE names
#    awk -v mapfile="${base}_peak_gene_class.tsv" '
#    BEGIN {
#        OFS = "\t"
#        while ((getline < mapfile) > 0) {
#            peak = $1
#            gene_id[peak] = $2
#            gene_sym[peak] = $3
#            region_class[peak] = $4
#        }
#    }
#    {
#        peak = $4
#        key = $1":"$2":"$3":"peak
#        peak_len = $3 - $2
#        te_name = ($8 != ".") ? $8 : ""
#        overlap_bp[key] += $NF
#        if (te_name != "") {
#            if (!seen[key, te_name]) {
#                seen[key, te_name] = 1
#                te_names[key] = (te_names[key] == "") ? te_name : te_names[key] ";" te_name
#            }
#        }
#        coords[key] = $1"\t"$2"\t"$3
#        pid[key] = peak
#    }
#    END {
#        print "Chr", "Start", "End", "GeneID", "GeneSymbol", "GeneRegion_Classification", "TE_Overlap_Percent", "TE_Overlap_Bin", "Overlapping_TE_Names"
#        for (k in coords) {
#            ov = overlap_bp[k] + 0
#            split(coords[k], a, "\t")
#            len = a[3] - a[2]
#            pct = (len > 0) ? sprintf("%.2f", (ov / len) * 100) : "0.00"
#            bin = (pct == 0) ? "0%" :
#                  (pct <= 10) ? "<=10%" :
#                  (pct <= 25) ? "<=25%" :
#                  (pct <= 50) ? "<=50%" :
#                  (pct <= 75) ? "<=75%" :
#                  (pct < 100) ? "<100%" : "100%"
#            gene = gene_id[pid[k]]
#            sym = gene_sym[pid[k]]
#            cls = region_class[pid[k]]
#            testr = (te_names[k] != "") ? te_names[k] : "None"
#            print coords[k], gene, sym, cls, pct, bin, testr
#        }
#    }' ${base}_TEraw.txt > "$OUT_DIR/${base}_TE_summary.tsv"

    # Cleanup temporary files
#    rm -f ${base}.bed ${base}_TEraw.txt ${base}_peak_gene_class.tsv

#done

#echo "Done. All summaries saved to: $OUT_DIR"

##6.6.25 trying to get gene symbols (Successful)
#module load BEDTools/2.31.0-GCC-12.3.0

# ==== SET PATHS ====
#ANN_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaksnew/ann"
#TE_BED="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/TEann_35_0.1filt.bed"
#GTF="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/refann.gtf"
#OUT_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaksnew/SummaryTables1"
#mkdir -p "$OUT_DIR"

# Create ENSDARG to GeneSymbol map from GTF
#awk -F'\t' '$3 == "gene" && $9 ~ /gene_id/ && $9 ~ /gene_name/ {
#    match($9, /gene_id "([^\"]+)"/, gid)
#    match($9, /gene_name "([^\"]+)"/, gname)
#    if (gid[1] && gname[1]) print gid[1] "\t" gname[1]
#}' "$GTF" > /tmp/ensid_to_symbol.tsv

# ==== PROCESS EACH PEAK FILE ====
#for annfile in $ANN_DIR/*.1000bp_ann.txt $ANN_DIR/*.5000bp_ann.txt; do
#    base=$(basename "$annfile" .txt)
#    echo "Processing $base"

    # Generate temp BED with PeakID
#    awk 'NR>1 {OFS="\t"; print $2, $3, $4, "peak"NR}' "$annfile" > ${base}.bed

    # Extract PeakID to GeneID and exon/intron classification
#    awk -F'\t' 'NR > 1 {
#        peak = "peak"NR
#        gene_id = $12
#        exon = ($9 ~ /exon/) ? 1 : 0
#        intron = ($9 ~ /intron/) ? 1 : 0
#        class = (exon && !intron) ? "exon_only" :
#                (!exon && intron) ? "intron_only" :
#                (exon && intron) ? "exon+intron" : "unclassified"
#        print peak, gene_id, class
#    }' "$annfile" > ${base}_peak_gene_class.tsv

    # Intersect with TE BED (include non-overlapping)
#    bedtools intersect -a ${base}.bed -b "$TE_BED" -wao > ${base}_TEraw.txt

    # Summarize TE overlap: percent, bin, and TE names with gene symbol mapping
#    awk -v mapfile="${base}_peak_gene_class.tsv" -v genesym="/tmp/ensid_to_symbol.tsv" '
#    BEGIN {
#        OFS = "\t"
#        while ((getline < mapfile) > 0) {
#            peak = $1
#            gene_id[peak] = $2
#            region_class[peak] = $3
#        }
#        while ((getline < genesym) > 0) {
#            id2sym[$1] = $2
#        }
#    }
#    {
#        peak = $4
#        key = $1":"$2":"$3":"peak
#        peak_len = $3 - $2
#        te_name = ($8 != ".") ? $8 : ""
#        overlap_bp[key] += $NF
#        if (te_name != "") {
#            if (!seen[key, te_name]) {
#                seen[key, te_name] = 1
#                te_names[key] = (te_names[key] == "") ? te_name : te_names[key] ";" te_name
#            }
#        }
#        coords[key] = $1"\t"$2"\t"$3
#        pid[key] = peak
#    }
#    END {
#        print "Chr", "Start", "End", "GeneID", "GeneSymbol", "GeneRegion_Classification", "TE_Overlap_Percent", "TE_Overlap_Bin", "Overlapping_TE_Names"
#        for (k in coords) {
#            ov = overlap_bp[k] + 0
#            split(coords[k], a, "\t")
#            len = a[3] - a[2]
#            pct = (len > 0) ? sprintf("%.2f", (ov / len) * 100) : "0.00"
#            bin = (pct == 0) ? "0%" :
#                  (pct <= 10) ? "<=10%" :
#                  (pct <= 25) ? "<=25%" :
#                  (pct <= 50) ? "<=50%" :
#                  (pct <= 75) ? "<=75%" :
#                  (pct < 100) ? "<100%" : "100%"
#            gid = gene_id[pid[k]]
#            sym = (gid in id2sym) ? id2sym[gid] : "NA"
#            cls = region_class[pid[k]]
#            testr = (te_names[k] != "") ? te_names[k] : "None"
#            print coords[k], gid, sym, cls, pct, bin, testr
#        }
#    }' ${base}_TEraw.txt > "$OUT_DIR/${base}_TE_summary.tsv"

    # Cleanup temporary files
#    rm -f ${base}.bed ${base}_TEraw.txt ${base}_peak_gene_class.tsv

#done

#rm /tmp/ensid_to_symbol.tsv

#echo "Done. All summaries saved to: $OUT_DIR"

##6.6.25 Trying to get 0%TE H3K9me3 peaks at genes. This did not work properly gave a list with a bunch of TEs 
# Directory where your TE summary tables are stored
#SUMMARY_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaksnew/SummaryTables1"

# Output file
#OUTPUT="${SUMMARY_DIR}/genes_with_H3K9me3_noTE_all.tsv"

# Header
#echo -e "GeneID\tGeneSymbol\tTimepoint_Window" > "$OUTPUT"

# Loop through each summary file
#for file in "$SUMMARY_DIR"/*_TE_summary.tsv; do
    # Extract just the filename, e.g., 3hpf_K9.5000bp_ann_TE_summary.tsv
#    fname=$(basename "$file")

    # Get source label like: 3hpf_K9_5000bp
#    label=$(echo "$fname" | sed -E 's/_TE_summary.tsv$//' | sed -E 's/\.ann/_ann/' | tr '.' '_')

    # Extract GeneID, GeneSymbol, and Source label for rows with 0% overlap
#    awk -v lbl="$label" -F'\t' 'NR > 1 && $7 == "0.00" { print $4 "\t" $5 "\t" lbl }' "$file"
#done | sort -u >> "$OUTPUT"

#echo "Done: Combined output saved to $OUTPUT"

##6.6.25 Trying to get files that will show the right genes with H3K9me3 enrichment. Homer seems to be misannotating introns, exons, etc. Since slitrk4 that we know is exon enriched, doesn't come up in this output
#module load BEDTools/2.31.0-GCC-12.3.0

#ANN_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaksnew/ann"
#TE_BED="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaks/TEann_35_0.1filt.bed"
#GTF="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/refann.gtf"
#OUT_DIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published/peaksnew/SummaryTables3"
#mkdir -p "$OUT_DIR"

# === Build gene symbol map from GTF ===
#awk -F'\t' '$3 == "gene" && $9 ~ /gene_id/ && $9 ~ /gene_name/ {
#    match($9, /gene_id "([^"]+)"/, gid)
#    match($9, /gene_name "([^"]+)"/, gname)
#    if (gid[1] && gname[1]) print gid[1] "\t" gname[1]
#}' "$GTF" > /tmp/ensid_to_symbol.tsv

# === Loop through peak annotation files ===
#for annfile in $ANN_DIR/*.1000bp_ann.txt $ANN_DIR/*.5000bp_ann.txt; do
#    base=$(basename "$annfile" .txt)
#    echo "Processing $base..."

    # Step 1: BED from HOMER peaks with peak ID
#    awk 'NR > 1 {OFS="\t"; print $2, $3, $4, "peak"NR}' "$annfile" > ${base}.bed

    # Step 2: extract PeakID → GeneID + Exon/Intron Class
#    awk -F'\t' '
#    NR > 1 {
#        peak = "peak"NR
#        gene_id = $12
#        exon = ($9 ~ /exon/) ? 1 : 0
#        intron = ($9 ~ /intron/) ? 1 : 0
#        class = (exon && !intron) ? "exon_only" :
#                (!exon && intron) ? "intron_only" :
#                (exon && intron) ? "exon+intron" : "unclassified"
#        print peak, gene_id, class
#    }' "$annfile" > ${base}_peak_to_gene.tsv

    # Step 3: bedtools intersect with TEs (with -wao)
#    bedtools intersect -a ${base}.bed -b "$TE_BED" -wao > ${base}_TEraw.txt

    # Step 4: Aggregate per GENE
#    awk -v pk2gene="${base}_peak_to_gene.tsv" -v symtab="/tmp/ensid_to_symbol.tsv" '
#    BEGIN {
#        OFS="\t"
#        while ((getline < pk2gene) > 0) {
#            peak = $1; gid = $2; class = $3
#            peak2gene[peak] = gid
#            peakclass[peak] = class
#        }
#        while ((getline < symtab) > 0) {
#            id2sym[$1] = $2
#        }
#    }
#    {
#        peak = $4
#        gid = peak2gene[peak]
#        if (gid == "") next
#        class = peakclass[peak]
#        geneclass[gid][class] = 1
#        genepeaklen[gid][peak] = $3 - $2
#        geneoverlap[gid][peak] += $NF
#        if ($8 != "." && $8 != "") {
#            if (!seen[gid, $8]) {
#                seen[gid, $8] = 1
#                tenames[gid] = (tenames[gid] == "") ? $8 : tenames[gid] ";" $8
#            }
#        }
#    }
#    END {
#        print "GeneID", "GeneSymbol", "GeneRegion_Classification", "TE_Overlap_Percent", "TE_Overlap_Bin", "Overlapping_TE_Names"
#        for (gid in genepeaklen) {
            # Aggregate classification
#            exon = intron = 0
#            for (cls in geneclass[gid]) {
#                if (cls ~ /exon/) exon = 1
#                if (cls ~ /intron/) intron = 1
#            }
#            region_class = (exon && intron) ? "exon+intron" :
#                           (exon && !intron) ? "exon_only" :
#                           (!exon && intron) ? "intron_only" : "unclassified"

#            # Max % overlap across all peaks
#            max_pct = 0
#            for (pk in genepeaklen[gid]) {
#                ov = geneoverlap[gid][pk] + 0
#                len = genepeaklen[gid][pk]
#                pct = (len > 0) ? (ov / len) * 100 : 0
#                if (pct > max_pct) max_pct = pct
#            }

            # Bin
#            bin = (max_pct == 0) ? "0%" :
#                  (max_pct <= 10) ? "<=10%" :
#                  (max_pct <= 25) ? "<=25%" :
#                  (max_pct <= 50) ? "<=50%" :
#                  (max_pct <= 75) ? "<=75%" :
#                  (max_pct < 100) ? "<100%" : "100%"

#            sym = (gid in id2sym) ? id2sym[gid] : "NA"
#            tename = (tenames[gid] != "") ? tenames[gid] : "None"

#            printf "%s\t%s\t%s\t%.2f\t%s\t%s\n", gid, sym, region_class, max_pct, bin, tename
#        }
#    }' ${base}_TEraw.txt > "$OUT_DIR/${base}_TE_byGene.tsv"

    # Cleanup
#    rm -f ${base}.bed ${base}_TEraw.txt ${base}_peak_to_gene.tsv
#done

#rm /tmp/ensid_to_symbol.tsv
#echo "All gene-level summaries saved to: $OUT_DIR"

##6.6.25 This is to reclassify the peaks with the gtf file that is good
module load BEDTools/2.31.0-GCC-12.3.0
BASEDIR="/scratch/dr27977/H3K9me3_Zebrafish/CUTnRUN_published"
GTF="${BASEDIR}/refann.gtf"
BEDDIR="${BASEDIR}/peaksnew/ann"
OUTDIR="${BASEDIR}/SummaryTables4"
CLASSIFIED_DIR="${BASEDIR}/GTF_Classified1"

mkdir -p "$OUTDIR"
mkdir -p "$CLASSIFIED_DIR"

# === Extract GTF annotations ===
awk '$3 == "gene"' "$GTF" | awk '{
    match($0, /gene_id "([^"]+)"/, a);
    match($0, /gene_name "([^"]+)"/, b);
    print $1"\t"$4-1"\t"$5"\t"a[1]"\t0\t"$7"\t"b[1]
}' > "${CLASSIFIED_DIR}/all_genes.bed"

awk '$3 == "exon"' "$GTF" | awk '{
    match($0, /gene_id "([^"]+)"/, a);
    print $1"\t"$4-1"\t"$5"\t"a[1]"\t0\t"$7
}' > "${CLASSIFIED_DIR}/all_exons.bed"

awk '$3 == "transcript"' "$GTF" | awk '{
    match($0, /gene_id "([^"]+)"/, a);
    print $1"\t"$4-1"\t"$5"\t"a[1]"\t0\t"$7
}' | bedtools sort -i - | bedtools subtract -a - -b "${CLASSIFIED_DIR}/all_exons.bed" > "${CLASSIFIED_DIR}/all_introns.bed"

# Sort all
sort -k1,1 -k2,2n "${CLASSIFIED_DIR}/all_genes.bed" > "${CLASSIFIED_DIR}/all_genes.sorted.bed"
sort -k1,1 -k2,2n "${CLASSIFIED_DIR}/all_exons.bed" > "${CLASSIFIED_DIR}/all_exons.sorted.bed"
sort -k1,1 -k2,2n "${CLASSIFIED_DIR}/all_introns.bed" > "${CLASSIFIED_DIR}/all_introns.sorted.bed"

# === Process all peak files ===
for peakbed in "${BEDDIR}"/*ann.bed; do
    base=$(basename "$peakbed" .bed)

    exon_hits="${CLASSIFIED_DIR}/${base}_has_exon.bed"
    intron_hits="${CLASSIFIED_DIR}/${base}_has_intron.bed"
    gene_both="${CLASSIFIED_DIR}/${base}_both.bed"
    gene_exon_only="${CLASSIFIED_DIR}/${base}_exon_only.bed"
    gene_intron_only="${CLASSIFIED_DIR}/${base}_intron_only.bed"

    # Intersect with exon and intron
    bedtools intersect -wa -a "$peakbed" -b "${CLASSIFIED_DIR}/all_exons.sorted.bed" > "$exon_hits"
    bedtools intersect -wa -a "$peakbed" -b "${CLASSIFIED_DIR}/all_introns.sorted.bed" > "$intron_hits"

    # Sort and uniq to avoid duplicated regions
    sort -u "$exon_hits" > "${exon_hits}.uniq"
    sort -u "$intron_hits" > "${intron_hits}.uniq"

    # Classify overlaps
    bedtools intersect -a "${exon_hits}.uniq" -b "${intron_hits}.uniq" > "$gene_both"
    bedtools intersect -v -a "${exon_hits}.uniq" -b "${intron_hits}.uniq" > "$gene_exon_only"
    bedtools intersect -v -a "${intron_hits}.uniq" -b "${exon_hits}.uniq" > "$gene_intron_only"

    # Add headers to table
    outfile="${OUTDIR}/${base}_peak_classification.tsv"
    echo -e "Gene_ID\tGene_Symbol\tClassification\tPeak_File" > "$outfile"

    for file in "$gene_exon_only" "$gene_intron_only" "$gene_both"; do
        class="unknown"
        case $file in
            *exon_only.bed) class="exon_only" ;;
            *intron_only.bed) class="intron_only" ;;
            *both.bed) class="both" ;;
        esac
        awk -v c="$class" -v pf="$base" '{print $4 "\t" $7 "\t" c "\t" pf}' "$file" >> "$outfile"
    done
done