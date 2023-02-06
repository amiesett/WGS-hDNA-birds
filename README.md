# WGS-hDNA-birds
Steps to analyzing sequencing and lab data from 60 shallow shotgun sequencing libraries. Libraries were prepared from 20 DNA extractions of 8 *Turdus* specimens collected in 1926/1927.

1-7 require basename_array.csv

1. Trim polynucleotide tail from read 2 of IDT libraries via seqtk
        01_trim_IDT.sh
2. Filter out PCR and optical duplicate via hts_SuperDeduper
        02_dedup.sh
3. Trim adapters and merge reads as necessary via SeqPrep; trim low quality bases from reads ends and filter short reads via Trimmomatic; remove low complexity reads using remove_low_complex.py from NF-polish pipeline
        03_cleanv2.sh
4. Align cleaned reads to the reference and index mapped reads via BWA and Samtools, respectively.
        04_align-merged.sh
5. Run mapDamage to assess influence of age on sequence data
        05_map_damage.sh
6. Collect stats from data filtering and trimming via seqkit; output read length distributions based on unclipped, mapped bases and other stats from mapped data via Samtools
        06_align_stats.sh
7. Output vcfs for each library based on differences from the reference genome to roughly assess contamination levels
        07_ref_difs.sh
8. Parse and analyze sequencing metrics from steps 1-7 and extraction and library data collected during labwork (including bioanalyzer traces). Requires:
   - bioanalyzer .xml files
   - Turdus_shallow_stats.csv - summary of sequencing metrics collected from 1-7
   - Turdus_shallow_rl_metrics.csv - summary of read length distributions
   - WGS-hDNA-birds-lab-data - extraction and library data collected during labwork

          WGS-hDNA-birds.R
