This folder contains all the script for human motif analysis and prep for TFBS analysis
human_motif_processing.R: process NME or MML data from human into each motif for DNase and non-DNase region for TFBS analysis
motif_break_array.R: function that running motifBreakR by splitting motifs into n segment and run each in
motifbreakR_parallel.R: modified version of motifBreakR that fix the issue of parallel in linux
motif_human_allelic.R: function that perform the binomial test for each motif on whether it's associated with NME or MML and annotate motif into known OMIM disease 
motif_array_all_JASPAR_default.sh: bash script to submit the arrary of jobs for motifbreakR using default settings
parsing_GOrilla.R: R script that reformat the table from GOrilla into data S5