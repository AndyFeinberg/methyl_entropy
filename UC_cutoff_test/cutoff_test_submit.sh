#!/bin/bash
#$ -l cee,mf=15G,h_vmem=20G
#$ -pe local 5
#$ -wd ../
#$ -o "../logfiles/cutoff-test-$JOB_ID.err"
#$ -j y
#$ -m ebas
#$ -M "yfang27@jhmi.edu"
#$ -N cutoff_test
# Args
cutoff="$1"
echo "${cutoff}"

# Command for parallelization within a node
Rscript --vanilla UC_cutoff_test/cutoff_test.R ${cutoff}