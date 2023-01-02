#!/bin/bash -l
#SBATCH --time=5:0:0
#SBATCH --mem=200G
#SBATCH --partition=lrgmem
ml R
Rscript make.R $1
