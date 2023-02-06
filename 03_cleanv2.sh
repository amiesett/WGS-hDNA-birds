#!/bin/bash
#SBATCH --job-name=03_cleanv2
#SBATCH --output=/lustre/scratch/asettlec/century_DNA/final/logs/%x_%a-%A.o
#SBATCH --mail-user=asettlec@ttu.edu
#SBATCH --mail-type=ALL
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=4
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=5300M
#SBATCH --array=1-64
date

source /home/asettlec/miniconda3/etc/profile.d/conda.sh
conda activate nf-polish

# define main working directory and threads to use
# *indir is raw fastqs when not including unpaired R1s*
indir=$SCRATCH/century_DNA/02_deduped
outdir=$SCRATCH/century_DNA/final/03_cleaned
array=$SCRATCH/century_DNA/basename_array.csv
minq=15

#########################################
# Should not need to edit past this point
#########################################

# Assign sample to this Slurm task via slurm array
basename=$( head -n${SLURM_ARRAY_TASK_ID} $array | tail -n1 )

##################################
# SeqPrep (adapters and merging)
##################################
cd ${outdir}

SeqPrep -f ${indir}/${basename}_dedup_R1.fastq.gz \
-r ${indir}/${basename}_dedup_R2.fastq.gz \
-1 ${basename}_adapt_merge_R1.fastq.gz \
-2 ${basename}_adapt_merge_R2.fastq.gz \
-A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
-3 ${basename}_discard_R1.fastq.gz \
-4 ${basename}_discard_R2.fastq.gz \
-s ${basename}_adapt_merge.fastq.gz

##################################
# Trimmomatic (quality)
##################################

# Trim paired reads based on quality via trimmomatic
trimmomatic PE -phred33 \
       -threads $SLURM_NTASKS \
       ${outdir}/${basename}_adapt_merge_R1.fastq.gz \
       ${outdir}/${basename}_adapt_merge_R2.fastq.gz \
       ${outdir}/${basename}_qual_R1.fastq.gz \
       ${outdir}/${basename}_qual_1u.fastq.gz \
       ${outdir}/${basename}_qual_R2.fastq.gz \
       ${outdir}/${basename}_qual_2u.fastq.gz \
       LEADING:$minq TRAILING:$minq SLIDINGWINDOW:4:$minq MINLEN:30

# Trim unpaired and merged reads based on quality via trimmomatic
trimmomatic SE -phred33 \
       -threads $SLURM_NTASKS \
       ${outdir}/${basename}_adapt_merge.fastq.gz \
       ${outdir}/${basename}_qual_U.fastq.gz \
       LEADING:$minq TRAILING:$minq SLIDINGWINDOW:4:$minq MINLEN:30

# Comment out unless including unpaired R1s
# Concatenate newly unpaired R1s to merged read file
cat ${outdir}/${basename}_qual_1u.fastq.gz >> ${outdir}/${basename}_qual_U.fastq.gz

##################################
# Low complexity removal
##################################

# Remove low complexity reads w/ nf-polish custom python script
# * -1 and -2 directories are "indir" when including unpaired R1s
remove_low_complex.py \
   -1 ${outdir}/${basename}_qual_R1.fastq.gz \
   -2 ${outdir}/${basename}_qual_R2.fastq.gz \
   -u ${outdir}/${basename}_qual_U.fastq.gz \
   -c 0.5 \
   -p ${basename}

 mv ${basename}_R1.fastq ${basename}_final_R1.fastq
 mv ${basename}_R2.fastq ${basename}_final_R2.fastq
 mv ${basename}_U.fastq ${basename}_final_U.fastq

# Compress data output as fastq to fastq.gz
pigz -p $SLURM_NTASKS -v -r ${basename}_final_R1.fastq
pigz -p $SLURM_NTASKS -v -r ${basename}_final_R2.fastq
pigz -p $SLURM_NTASKS -v -r ${basename}_final_U.fastq

date
