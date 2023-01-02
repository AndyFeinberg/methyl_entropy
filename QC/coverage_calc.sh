#!/bin/bash
#SBATCH -J coverage_calc
#SBATCH --partition=shared,parallel,skylake
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH -o coverage_logfile/coverage-%A.out
#SBATCH -e coverage_logfile/coverage-%A.err
bed_in=$1
bam_in=$2
cov_out=$3


bedtools coverage -f 1 -sorted -counts -a ${bed_in} -b ${bam_in} > ${cov_out}
