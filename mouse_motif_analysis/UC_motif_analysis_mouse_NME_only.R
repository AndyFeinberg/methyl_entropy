library(data.table)
library(dplyr)
library(glue)
library(GenomicRanges)
library(parallel)

##function for calculating known motif enrichment
enrich_motif <- function(
	target_region,
	control_region,
	functional_region = NULL,
	region_size = 200,
	motif_lib = "JASPAR_hg19"
){
	if(class(target_region)[1] != "GRanges"){
		colnames(target_region)[1:3] <- c("chr","start","end")
		target_region <- makeGRangesFromDataFrame(target_region)
	}

	if(class(control_region)[1] != "GRanges"){
		colnames(control_region)[1:3] <- c("chr","start","end")
		control_region <- makeGRangesFromDataFrame(control_region)
	}

	target_region <- resize(target_region, width = region_size, fix='center')
	control_region <- resize(control_region, width = region_size, fix='center')

	switch(motif_lib,
		JASPAR_hg19 = {motif_lib_data <- motif_lib_hg19
						message("Using reference genome hg19")},
		JASPAR_hg38 = {motif_lib_data <- motif_lib_hg38
						message("Using reference genome hg38")},
		JASPAR_mm10 = {motif_lib_data <- motif_lib_mm10
						message("Using reference genome mm10")},
		stop("Please enter the correct reference genome")
	)

	if(!is.null(functional_region)){

		if(class(functional_region)[1] != "GRanges"){
			colnames(functional_region)[1:3] <- c("chr","start","end")
			functional_region <- makeGRangesFromDataFrame(functional_region)
		}
		message(glue("Filtering motif site using {length(functional_region)} functional regions"))

		motif_lib_data <- lapply(motif_lib_data, function(x){
			overlap_motif <- findOverlaps(x, functional_region) %>% queryHits()
			x[overlap_motif]
		})

		motif_lib_data <- GRangesList(motif_lib_data)
	}

	target_overlap <- countOverlaps(motif_lib_data, target_region)
	control_overlap <- countOverlaps(motif_lib_data, control_region)

	motif_hit <- data.frame(target_hit = target_overlap, control_hit = control_overlap)
	motif_off <- data.frame(target_off = length(target_region) - motif_hit$target_hit, control_off = length(control_region) - motif_hit$control_hit)
	row.names(motif_off) <- row.names(motif_hit)

	test_result <- lapply(c(1:nrow(motif_hit)), function(x){

		test_table <- matrix(c(motif_hit[x,1]+1,motif_hit[x,2]+1,motif_off[x,1]+1,motif_off[x,2]+1),
	                         nrow = 2,dimnames = list(c("target", "control"),c("hit", "off")))
		odds_se <- sqrt(sum(1/test_table))
	  test_result <- fisher.test(test_table, alternative = "greater")
		pval <- test_result$p.value
		odds <- test_result$estimate
		return(data.frame(odds_ratio=odds, pvalue=pval, odds_ratio_se=odds_se))
	})

	test_result <- Reduce(rbind, test_result)
	row.names(test_result) <- row.names(motif_hit)

	test_result$FDR <- p.adjust(test_result$pvalue, method="BH")

	return(cbind(test_result, motif_hit))
}

##load motif site database
motif_lib_mm10 <- readRDS("../downstream/input/mouse_analysis/motif_analysis/motif_JASPAR_mm10.rds")

setwd("../downstream/output/mouse_analysis/UC_motif_enrichment")

##sub region type: "NME only", "Both", "Neither", "MML only"
cluster_gr_anno <- readRDS(file="../downstream/input/mouse_analysis/tissue_region_enhancer_gr.rds")
cluster_gr_control <- readRDS(file="../downstream/input/mouse_analysis/tissue_region_control_enhancer_gr.rds")

tissue_name <- names(cluster_gr_anno)

##calculate motif enrichment for NME only regions
motif_enrichment_tissue <- mclapply(tissue_name, function(x){

	data_test <- cluster_gr_anno[[x]]
	target <- data_test[data_test$region_type == "NME only"]
	control <- cluster_gr_control[[x]]

	motif_enrichment <- enrich_motif(
		target_region = target,
		control_region = control,
		region_size = 250,
		motif_lib = "JASPAR_mm10"
	)

	saveRDS(motif_enrichment, file=glue("./enhancer_region/{x}_motif_enrichment_DNase_control_dNME.rds"))
},mc.cores=4)

q(save="no")
