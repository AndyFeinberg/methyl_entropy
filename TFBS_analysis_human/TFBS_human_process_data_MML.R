library(GenomicRanges)
library(dplyr)
library(glue)
library(matrixStats)

setwd("../")

##read in MML data for regulatory DNA regions and calculate the median MML for each motif
data_median_DNase <- lapply(c(1:3), function(x) {

	data_temp <- readRDS(glue("../downstream/output/human_analysis/Ken_motif/homogeneous/JASPAR_motif_hg19_MML_{x}_agnostic_DNase.rds"))
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
saveRDS(data_median_com_DNase, file="../downstream/output/human_analysis/Ken_motif/homogeneous/human_motif_DNase_median_MML.rds")

##read in MML data for non-regulatory DNA regions and calcualte the median MML for each motif
data_median_control <- lapply(c(1:3), function(x) {

	data_temp <- readRDS(glue("../downstream/output/human_analysis/Ken_motif/homogeneous/JASPAR_motif_hg19_MML_{x}_agnostic_Control.rds"))
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
saveRDS(data_median_com_control, file="../downstream/output/human_analysis/Ken_motif/homogeneous/human_motif_control_median_MML.rds")

q(save="no")
