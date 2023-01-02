library(GenomicRanges)
library(dplyr)
library(glue)
library(matrixStats)

setwd("../")

##read in NME data for regulatory DNA regions and calculate the median NME for each motif
data_median_DNase <- mclapply(c(1:3), function(x) {

	data_temp <- readRDS(glue("../downstream/output/human_analysis/Ken_motif/homogeneous/JASPAR_motif_hg19_NME_{x}_agnostic_DNase.rds"))
	data_temp_median <- sapply(data_temp, function(y){
		data_mat <- as.matrix(mcols(y))
		data_mat_median <- colMedians(data_mat,na.rm = TRUE)
		names(data_mat_median) <- colnames(data_mat)
		return(data_mat_median)
	})
	data_temp_median <- t(data_temp_median)

	return(data_temp_median)
},mc.cores=3)

data_median_com_DNase <- Reduce(rbind, data_median_DNase)
saveRDS(data_median_com_DNase, file=hg19_DNase_NME_fn)

##read in NME data for non-regulatory DNA regions and calcualte the median NME for each motif
data_median_control <- mclapply(c(1:3), function(x) {

	data_temp <- readRDS(glue("../downstream/output/human_analysis/Ken_motif/homogeneous/JASPAR_motif_hg19_NME_{x}_agnostic_Control.rds"))
	data_temp_median <- sapply(data_temp, function(y){
		data_mat <- as.matrix(mcols(y))
		data_mat_median <- colMedians(data_mat,na.rm = TRUE)
		names(data_mat_median) <- colnames(data_mat)
		return(data_mat_median)
	})
	data_temp_median <- t(data_temp_median)

	return(data_temp_median)
},mc.cores=3)

data_median_com_control <- Reduce(rbind, data_median_control)
saveRDS(data_median_com_control, file=hg19_control_NME_fn)

q(save="no")
