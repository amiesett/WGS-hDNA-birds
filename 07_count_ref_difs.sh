#!/bin/bash
#SBATCH --job-name=07_count_ref_difs
#SBATCH --output=/lustre/scratch/asettlec/century_DNA/final/logs/%x_%a-%A.o
#SBATCH --mail-user=asettlec@ttu.edu
#SBATCH --mail-type=ALL
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=5300M
#SBATCH --array=1-64
date

module load intel/18.0.3.222 samtools/1.9 bcftools/1.9

# define main working directory and threads to use
indir=$SCRATCH/century_DNA/final/04_bams
outdir=$SCRATCH/century_DNA/final/07_vcfs
array=$SCRATCH/century_DNA/basename_array.csv
reference=/home/jmanthey/references/GCA_013186435.1_ASM1318643v1_genomic.fna

#########################################
# Should not need to edit past this point
#########################################

basename=$( head -n${SLURM_ARRAY_TASK_ID} $array | tail -n1 )

# mpileup vcf
bcftools mpileup -Ou -f $reference ${indir}/${basename}_final.bam | bcftools call -vmO v -o ${outdir}/${basename}.vcf

# textfile of 1 allele call per line
grep -v ^# ${outdir}/${basename}.vcf | cut -f 10 | cut -d":" -f1 | sed 's/\//\n/' > ${outdir}/${basename}_gt_calls.txt

date
