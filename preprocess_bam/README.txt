The code in this folder is to process downloaded bam file from mouse
deduplicate.sh: remove duplicate for target bam file
merge_bam.sh: merge the bam file between replicates
soft_clip_bam.sh: remove 5bp from begining and end of each read like CPEL did for plotting, the clipped file is NOT send to CPEL, only the unclipped ones are sent to CPEL.