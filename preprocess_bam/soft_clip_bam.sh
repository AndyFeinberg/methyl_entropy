#!/bin/bash
#SBATCH -J trimbam
#SBATCH --partition=shared,skylake,express,parallel
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=10
#SBATCH -o trim_logfile/trim-%A.out
#SBATCH -e trim_logfile/trim-%A.err
bam_in=$1
bam_out=/ibox/afeinbe2/yfang/trim_bam/$(echo $bam_in|sed 's/\.bam//')
echo processing: $bam_in
echo output name: $bam_out
bam trimbam $bam_in $bam_out.trimmed.sam -L 5 -R 5 -c
samtools view -@10 -b -o $bam_out.trimmed.bam $bam_out.trimmed.sam
rm $bam_out.trimmed.sam
samtools index -@10 $bam_out.trimmed.bam
