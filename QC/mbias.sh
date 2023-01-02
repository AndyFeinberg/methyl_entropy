#!/bin/bash
#SBATCH -J mbias
#SBATCH --partition=shared,parallel,skylake
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH -o mbias_logfile/mbias-%A.out
#SBATCH -e mbias_logfile/mbias-%A.err
bam_in=$1
out_dir=$2
genome_folder=$3
echo processing: $bam_in
bismark_methylation_extractor -o $out_dir --parallel 12 --genome_folder $genome_folder --mbias_only -s $bam_in

