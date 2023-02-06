#!/bin/bash
#SBATCH --job-name=06_align_stats
#SBATCH --output=/lustre/scratch/asettlec/century_DNA/final/logs/%x_%a-%A.o
#SBATCH --mail-user=asettlec@ttu.edu
#SBATCH --mail-type=ALL
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=5300M
#SBATCH --array=1-64
date

source /home/asettlec/miniconda3/etc/profile.d/conda.sh
conda activate nf-polish
module load intel/18.0.3.222 java/1.8.0 bwa/0.7.17 samtools/1.9

# define main working directory and threads to use
indir=$SCRATCH/century_DNA/final
outdir=$SCRATCH/century_DNA/final/06_stats
array=$SCRATCH/century_DNA/basename_array.csv

#########################################
# Should not need to edit past this point
#########################################

basename=$( head -n${SLURM_ARRAY_TASK_ID} $array | tail -n1 )

# Output read length distributions via Samtools and awk
# Mapped, paired R1s
samtools view -f 65 -F 4 ${indir}/04_bams/${basename}_final.bam | awk '{print length($10)}' | sort -h | uniq -c > ${outdir}/${basename}_lengths1.csv
# Mapped, paired R2s
samtools view -f 129 -F 4 ${indir}/04_bams/${basename}_final.bam | awk '{print length($10)}' | sort -h | uniq -c > ${outdir}/${basename}_lengths2.csv
# Mapped, merged reads
samtools view -F 5 ${indir}/04_bams/${basename}_final.bam | awk '{print length($10)}' | sort -h | uniq -c > ${outdir}/${basename}_lengthsU.csv

# Collect summary stats of raw, cleaned, merged reads and bases
# raw reads read into hts_SuperDeduper
seqkit stats ${indir}/00_fastq/${basename}_R1.fastq.gz -T > ${outdir}/${basename}.temp
raw_reads=$(( $( grep "R1" ${outdir}/${basename}.temp | cut -f4 ) * 2 ))
# raw base read into hts_SuperDeduper
raw_bases=$(( $( grep "R1" ${outdir}/${basename}.temp | cut -f5 ) * 2 ))

# reads output from hts_SuperDeduper
seqkit stats ${indir}/00_fastq/${basename}_R1.fastq.gz -T > ${outdir}/${basename}.temp
duplicate_reads=$(( $raw_reads - ($( grep "R1" ${outdir}/${basename}.temp | cut -f4 ) * 2) ))
# bases output from hts_SuperDeduper
duplicate_bases=$(( $raw_bases - ($( grep "R1" ${outdir}/${basename}.temp | cut -f5 ) * 2) ))

# cleaned reads
seqkit stats ${indir}/03_cleaned/${basename}_final_R1.fastq.gz ${indir}/03_cleaned/${basename}_final_R2.fastq.gz ${indir}/03_cleaned/${basename}_final_U.fastq.gz -T > ${outdir}/${basename}.temp
cleaned_reads=$(( $( grep "R1" ${outdir}/${basename}.temp | cut -f4 ) * 2 + $( grep "_U" ${outdir}/${basename}.temp | cut -f4 ) ))
# cleaned bases
cleaned_bases=$(( $( grep "R1" ${outdir}/${basename}.temp | cut -f5 ) + $( grep "R2" ${outdir}/${basename}.temp | cut -f5 ) + $( grep "_U" ${outdir}/${basename}.temp | cut -f5 ) ))

# Run samtools stats
samtools stats ${indir}/04_bams/${basename}_final.bam | grep ^SN | cut -f 2- > ${outdir}/${basename}.temp
# Collect summary stats from samtools stats output
mapped_reads=$(grep 'reads mapped:' ${outdir}/${basename}.temp | cut -f 2)
mapped_bases=$(grep 'bases mapped:' ${outdir}/${basename}.temp | cut -f 2)
mapped_unclipped_bases=$(grep 'bases mapped (cigar):' ${outdir}/${basename}.temp | cut -f 2)
mapped_paired=$(grep 'reads mapped and paired' ${outdir}/${basename}.temp | cut -f 2)
properly_paired=$(grep 'reads properly paired' ${outdir}/${basename}.temp | cut -f 2)
length=$(grep 'average length' ${outdir}/${basename}.temp | cut -f 2)
quality=$(grep 'average quality' ${outdir}/${basename}.temp | cut -f 2)
# Samtools stats def'n of insert size average - the average absolute template length for paired and mapped reads.
insert=$(grep 'insert size average' ${outdir}/${basename}.temp | cut -f 2)
insert_stdev=$(grep 'insert size standard deviation' ${outdir}/${basename}.temp | cut -f 2)
# Calculate average depth per basepair
depth=$(samtools depth ${indir}/04_bams/${basename}_final.bam |  awk '{sum+=$3} END { print sum }')
avg_depth=$(samtools depth ${indir}/04_bams/${basename}_final.bam |  awk '{sum+=$3} END { print sum/NR}')

# Output summary stats to file
echo -e "$basename,$raw_reads,$cleaned_reads,$mapped_reads,$raw_bases,$cleaned_bases,$mapped_bases,$mapped_unclipped_bases,$mapped_paired,$properly_paired,$length,$quality,$depth,$avg_depth" > ${outdir}/${basename}_stats.csv

rm ${outdir}/${basename}.temp

date

