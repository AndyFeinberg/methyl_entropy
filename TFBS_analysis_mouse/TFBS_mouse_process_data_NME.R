library(GenomicRanges)
library(dplyr)
library(glue)
library(matrixStats)

setwd("../downstream/output/mouse_analysis/TFBS_analysis")

##read in NME data for regulatory DNA regions and calculate the median NME for each motif
data_median_DNase <- lapply(c(1:12), function(x) {

	data_temp <- readRDS(glue("../downstream/input/mouse_analysis/TFBS_analysis/JASPAR_motif_mm10_NME_{x}_agnostic_merged_DNase.rds"))
	data_temp_median <- sapply(data_temp, function(y){
		data_mat <- as.matrix(mcols(y))
		data_mat_median <- colMedians(data_mat,na.rm = TRUE)
		names(data_mat_median) <- colnames(data_mat)
		return(data_mat_median)
	})
	data_temp_median <- t(data_temp_median)

	return(data_temp_median)
})

data_median_com_DNase <- Reduce(rbind, data_median_DNase)
saveRDS(data_median_com_DNase, file="mouse_motif_DNase_median.rds")

##read in NME data for non-regulatory DNA regions and calcualte the median NME for each motif
data_median_control <- lapply(c(1:12), function(x) {

	data_temp <- readRDS(glue("../downstream/input/mouse_analysis/TFBS_analysis/JASPAR_motif_mm10_NME_{x}_agnostic_merged_control.rds"))
	data_temp_median <- sapply(data_temp, function(y){
		data_mat <- as.matrix(mcols(y))
		data_mat_median <- colMedians(data_mat,na.rm = TRUE)
		names(data_mat_median) <- colnames(data_mat)
		return(data_mat_median)
	})
	data_temp_median <- t(data_temp_median)

	return(data_temp_median)
})

data_median_com_control <- Reduce(rbind, data_median_control)
saveRDS(data_median_com_control, file="mouse_motif_control_median.rds")

q(save="no")
