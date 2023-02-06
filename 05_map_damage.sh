#!/bin/bash
#SBATCH --job-name=05_map_damage_u1s
#SBATCH --output=/lustre/scratch/asettlec/century_DNA/v2-trim15q20/logs/%x_%a-%A.o
#SBATCH --mail-user=asettlec@ttu.edu
#SBATCH --mail-type=ALL
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=5300M
#SBATCH --array=5-64
date

source /home/asettlec/miniconda3/etc/profile.d/conda.sh
conda activate mapdamage2

# define main working directory and threads to use
indir=$SCRATCH/century_DNA/v2-trim15q20/04_bams
outdir=$SCRATCH/century_DNA/v2-trim15q20/05_damage
array=$SCRATCH/century_DNA/basename_array.csv
reference=/home/jmanthey/references/GCA_013186435.1_ASM1318643v1_genomic.fna

#########################################
# Should not need to edit past this point
#########################################

# Assign sample to this Slurm task via slurm array
basename=$( head -n${SLURM_ARRAY_TASK_ID} $array | tail -n1 )

##################################
# mapDamage
##################################
mkdir ${outdir}/temp_${basename}

mapDamage -i ${indir}/${basename}_final.bam -r $reference -d ${outdir}/temp_${basename} --merge-reference-sequences

cd ${outdir}/temp_${basename}
for f in $(ls); do mv $f ../${basename}_${f}; done
cd ../
rm -r temp_${basename}

date
