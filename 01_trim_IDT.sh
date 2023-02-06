#!/bin/bash
#SBATCH --job-name=01_trim_IDT10
#SBATCH --output=logs/%x_%a-%A.o
#SBATCH --mail-user=asettlec@ttu.edu
#SBATCH --mail-type=ALL
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=5300M
#SBATCH --array=1-64
date

source /home/asettlec/miniconda3/etc/profile.d/conda.sh
conda activate mapdamage2

# define main working directory and threads to use
workdir=$SCRATCH/century_DNA
indir=$SCRATCH/century_DNA/00_fastq/raw_IDT_R2s
outdir=$SCRATCH/century_DNA/00_fastq
array=$SCRATCH/century_DNA/basename_array.csv
trimn=10

#########################################
# Should not need to edit past this point
#########################################

# Assign sample to this Slurm task via slurm array
basename=$( head -n${SLURM_ARRAY_TASK_ID} $array | tail -n1 )

# Trim first 10 bp of R2 from each IDT library to remove adaptase tail
if [[ ${basename} = *IDT ]] || [[ ${basename} = *negative ]]
then
  seqtk trimfq -b $trimn ${indir}/${basename}_R2.fastq.gz > \
       ${outdir}/${basename}_trim${trimn}_R2.fastq
fi

date
