#!/bin/bash
#SBATCH -J deduplicate
#SBATCH --partition=shared,parallel,skylake
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=24
#SBATCH -o dedup_logfile/dedup-%A.out
#SBATCH -e dedup_logfile/dedup-%A.err
bam_in=$1
bam_in=$(echo $bam_in|sed 's/\.sort\.bam//')
echo processing: $bam_in
#(deduplicate_bismark -s --bam $bam_in.sort.bam)&& (samtools index $bam_in.sort.deduplicated.bam)
(java -Xmx95g -jar ~/data/yfang/bin/picard/build/libs/picard.jar MarkDuplicates INPUT=$bam_in.sort.bam OUTPUT=$bam_in.sort.dup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT QUIET=true METRICS_FILE=$bam_in.sort.dup_report.txt TMP_DIR=$PWD/picard_tmp READ_NAME_REGEX=null)&&(samtools index -@24 $bam_in.sort.dup.bam)
