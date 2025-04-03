#!/bin/bash
#SBATCH --job-name=trim		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1		                            # Single task job
#SBATCH --cpus-per-task=24		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=120gb			                            # Total memory for job
#SBATCH --time=48:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/dr27977/H3K9me3_ZF			# Location of standard output and error log files 
#SBATCH --mail-user=dr27977@uga.edu                    # Where to send mail 
#SBATCH --mail-type=BEGIN,END,FAIL,ALL                            # Mail events (BEGIN, END, FAIL, ALL)

OUTDIR="/scratch/dr27977/H3K9me3_ZF" 
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

for file in $BASEDIR/*_R*.fastq.gz;
do
  if [[ $prefix ]]; then
        base=$(basename ${first} _R1.fastq.gz)
        sh $HOMEDIR/PE_trim_and_star.sh -o $BASEDIR -n $base -m one $first $file
        prefix=
    else
        first=$file
        prefix=${file%%_*}
    fi
done

#Need to validate input files. Job script calls PE_trim_and_star, and Script is reading through files but is not making files 