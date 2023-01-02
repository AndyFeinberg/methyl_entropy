#!/bin/bash
#SBATCH -J merge_bam
#SBATCH --partition=shared,parallel,skylake
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=24
#SBATCH -o "merge_bam_log/merge_bam-%j.out"
#SBATCH -e "merge_bam_log/merge_bam-%j.err"
# Args
bam_in="$1"
echo samtools merge -@24 "$bam_in"_all.sort.dup.bam "$bam_in"_merged1.sort.dup.bam "$bam_in"_merged2.sort.dup.bam
echo samtools index -@24 "$bam_in"_all.sort.dup.bam
samtools merge -@24 "$bam_in"_all.sort.dup.bam "$bam_in"_merged1.sort.dup.bam "$bam_in"_merged2.sort.dup.bam
samtools index -@24 "$bam_in"_all.sort.dup.bam


