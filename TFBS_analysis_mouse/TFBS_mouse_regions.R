library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(exomeCopy)
options(scipen=999)

setwd("../downstream/input/mouse_analysis")

##obtain regulaotry DNA regions from mouse DNase-seq peak data
metadata <- readRDS("../downstream/input/mouse_analysis/ENCODE_data/metadata/match_RNA_DNase_mm10.rds")
metadata_filter <- metadata[grep("C57BL/6",metadata$Biosample),]
colnames(metadata_filter) <- c("sample", "RNA_ID", "DNase_ID")

match_file <- readRDS("../downstream/input/mouse_analysis/ENCODE_data/metadata/DNase_mm10_file_match.rds")
colnames(match_file) <- c("DNase_peak", "DNase_ID")
peak_match <- match_file[match_file$DNase_ID %in% metadata_filter$DNase_ID,]

peak_match <- left_join(peak_match, metadata_filter)
saveRDS(peak_match, file="DNase_RNA_peak_match_mm10.rds")

filelist <- list.files(path="../downstream/input/mouse_analysis/ENCODE_data/peak/mm10/DNase_house",pattern="\\.narrowPeak",full.names = TRUE, recursive = TRUE)

filelist_in <- sapply(peak_match[,1],function(x){
  grep(x,filelist, value = TRUE)
})

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
grl_narrowPeak <- lapply(filelist_in,function(x){
  import(x,format = "BED", extraCols = extraCols_narrowPeak)
})

grl_narrowPeak <- GRangesList(grl_narrowPeak)
grl_merge <- Reduce(c,grl_narrowPeak)

##merge peaks from all samples
grl_merge <- reduce(grl_merge)

select_chr <- c(paste0("chr",1:19),"chrX", "chrY")
grl_merge_filter <- grl_merge[seqnames(grl_merge) %in% select_chr]

##divide peak regions into 250bp bins
grl_merge_filter_bin <- subdivideGRanges(grl_merge_filter,subsize=250)

saveRDS(grl_merge_filter_bin,file="DNase_mm10_peak_merge_250bp.rds")

##run cisgenome to get control regions (i.e., non-regulatory DNA regions)
gr_peak_fix_cod <- data.frame(grl_merge_filter_bin)[,-4]
gr_peak_fix_cod$strand <- "+"

write.table(gr_peak_fix_cod, file="DNase_mm10_peak_merge_250bp.cod", col.names=FALSE, sep="\t", quote=FALSE)

system("~/cisgenome_project/bin/refgene_getmatchedcontrol -d ~/cisgenome_project/mm10/annotation/refFlat_sorted.txt -dt 1 -s mouse -c ~/cisgenome_project/mm10/chrlen.txt -i DNase_mm10_peak_merge_250bp.cod -o DNase_mm10_peak_merge_250bp_control.cod -n 1 -l 250 -nr 1")

##exclude overlapping regions from the control
peak_control <- fread("DNase_mm10_peak_merge_250bp_control.cod")
gr_peak_control <- makeGRangesFromDataFrame(peak_control[,c(2:4)])
gr_peak_control <- resize(gr_peak_control, 250, fix = "center")

gr_peak <- grl_merge_filter_bin

overlap <- findOverlaps(gr_peak_control, gr_peak)
gr_peak_control_filter <- gr_peak_control[-unique(queryHits(overlap))]

saveRDS(gr_peak_control_filter,file="DNase_mm10_peak_merge_250bp_control.rds")
