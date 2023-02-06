#!/bin/bash
#SBATCH --job-name=04_align-merged
#SBATCH --output=/lustre/scratch/asettlec/century_DNA/final/logs/%x_%a-%A.o
#SBATCH --mail-user=asettlec@ttu.edu
#SBATCH --mail-type=ALL
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=5300M
#SBATCH --array=1-64
date

module load intel/18.0.3.222 java/1.8.0 bwa/0.7.17 samtools/1.9

# define main working directory
indir=$SCRATCH/century_DNA/final/03_cleaned
outdir=$SCRATCH/century_DNA/final/04_bams
reference=/home/jmanthey/references/GCA_013186435.1_ASM1318643v1_genomic.fna
array=$SCRATCH/century_DNA/basename_array.csv

#########################################
# Should not need to edit past this point
#########################################

# Assign sample to this Slurm task via slurm array
basename=$( head -n${SLURM_ARRAY_TASK_ID} $array | tail -n1 )
# Limit threads and memory for some GATK processes that spinoff extra threads and memory to start
spark_threads=$(echo " ${SLURM_NTASKS} - 1" | bc -l)
spark_mem=$(echo " ${SLURM_MEM_PER_CPU} * 0.8" | bc -l)

# run bwa mem for unmerged reads
bwa mem -t $SLURM_NTASKS ${reference} ${indir}/${basename}_final_R1.fastq.gz ${indir}/${basename}_final_R2.fastq.gz > ${outdir}/${basename}_unmerged.sam

# run bwa mem for merged reads
bwa mem -t $SLURM_NTASKS ${reference} ${indir}/${basename}_final_U.fastq.gz > ${outdir}/${basename}_merged.sam

# convert sam to bam for unmerged reads
samtools view -b -S -o ${outdir}/${basename}_unmerged.bam ${outdir}/${basename}_unmerged.sam

# convert sam to bam for merged reads
samtools view -b -S -o ${outdir}/${basename}_merged.bam ${outdir}/${basename}_merged.sam

# combine the two bam files
samtools cat ${outdir}/${basename}_merged.bam ${outdir}/${basename}_unmerged.bam > ${outdir}/${basename}.bam

# remove ungrouped sam and bam files
rm ${outdir}/${basename}_unmerged.sam
rm ${outdir}/${basename}_unmerged.bam
rm ${outdir}/${basename}_merged.sam
rm ${outdir}/${basename}_merged.bam

# clean up the bam file
gatk CleanSam -I ${outdir}/${basename}.bam -O ${outdir}/${basename}_cleaned.bam

# remove the raw bam
rm ${outdir}/${basename}.bam

# add read groups to the cleaned bam file
gatk AddOrReplaceReadGroups -I ${outdir}/${basename}_cleaned.bam -O ${outdir}/${basename}_cleaned_rg.bam --RGLB 1 --RGPL illumina --RGPU unit1 --RGSM ${basename}

# remove the cleaned bam file
rm ${outdir}/${basename}_cleaned.bam

# Sort the cleaned and grouped bam file
gatk SortSam -I ${outdir}/${basename}_cleaned_rg.bam -O ${outdir}/${basename}_final.bam --SORT_ORDER coordinate

# Remove the unsorted bam file
rm ${outdir}/${basename}_cleaned_rg.bam

# index the final bam file -- edited command from JDM script to include -@ option to enable multithreading
samtools index -@ $SLURM_NTASKS ${outdir}/${basename}_final.bam

date
