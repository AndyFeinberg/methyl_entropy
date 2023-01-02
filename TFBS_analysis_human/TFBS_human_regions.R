library(GenomicRanges)
library(data.table)
options(scipen=999)

setwd("../downstream/input/human_analysis")

##read in DNase regions (i.e., regulatory DNA regions)
##DNase regions were obtained from https://github.com/WeiqiangZhou/BIRD-data/releases/download/v3.0/BIRD_data_ENCODE.zip and resize to 250bp for each region.
DNase_region <- readRDS("../downstream/input/human_analysis/DNase_hg19_250bp.rds")

##run cisgenome to get control regions (i.e., non-regulatory DNA regions)
gr_peak_fix_cod <- data.frame(DNase_region)[,-4]
gr_peak_fix_cod$strand <- "+"

write.table(gr_peak_fix_cod, file="DNase_hg19_250bp.cod", col.names=FALSE, sep="\t", quote=FALSE)

system("~/cisgenome_project/bin/refgene_getmatchedcontrol -d ~/cisgenome_project/hg19/annotation/refFlat_sorted.txt -dt 1 -s human -c ~/cisgenome_project/hg19/chrlen.txt -i DNase_hg19_250bp.cod -o DNase_hg19_250bp_control.cod -n 1 -l 250 -nr 1")

##exclude overlapping regions from the control
peak_control <- fread("DNase_hg19_250bp_control.cod")[-1,]
gr_peak_control <- makeGRangesFromDataFrame(peak_control[,c(2:4)])
gr_peak_control <- resize(gr_peak_control, 250, fix = "center")

gr_peak <- DNase_region

overlap <- findOverlaps(gr_peak_control, gr_peak)
gr_peak_control_filter <- gr_peak_control[-unique(queryHits(overlap))]

saveRDS(gr_peak_control_filter,file="DNase_hg19_250bp_control.rds")
