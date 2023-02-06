#!/bin/bash
#SBATCH --job-name=02_dedup_trim10
#SBATCH --output=logs/%x_%a-%A.o
#SBATCH --mail-user=asettlec@ttu.edu
#SBATCH --mail-type=ALL
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=5300M
#SBATCH --array=1-64
date

source /home/asettlec/miniconda3/etc/profile.d/conda.sh
conda activate nf-polish

# define main working directory and threads to use
indir=$SCRATCH/century_DNA/00_fastq
outdir=$SCRATCH/century_DNA/02_deduped
array=$SCRATCH/century_DNA/basename_array.csv
idt_trimn=10

#########################################
# Should not need to edit past this point
#########################################

# Assign sample to this Slurm task via slurm array
basename=$( head -n${SLURM_ARRAY_TASK_ID} $array | tail -n1 )

##################################
# SuperDeduper
##################################
mkdir ${outdir}/${basename}
cd ${outdir}/${basename}

# Run SuperDeduper to identify and remove duplicate reads
if [[ ${basename} = *IDT ]] || [[ ${basename} = *negative ]]
then
  hts_SuperDeduper -1 ${indir}/${basename}_R1.fastq.gz \
  -2 ${indir}/${basename}_trim${idt_trimn}_R2.fastq.gz -f PE
else
  hts_SuperDeduper -1 ${indir}/${basename}_R1.fastq.gz \
  -2 ${indir}/${basename}_R2.fastq.gz -f PE
fi

mv PE_R1.fastq.gz ../${basename}_dedup_R1.fastq.gz
mv PE_R2.fastq.gz ../${basename}_dedup_R2.fastq.gz
mv stats.log ../${basename}.stats
cd ..
rm -r ${outdir}/${basename}

date
