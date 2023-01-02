library(GenomicRanges)
library(dplyr)
library(glue)
library(ComplexHeatmap)
library(ggplot2)
library(patchwork)
library(data.table)

setwd("../downstream/output/mouse_analysis/UC_motif_enrichment")

region_type <- c("dMML", "dNME")
input_dir <- "enhancer_region"

##read in log odds ratio
result_type <- lapply(region_type, function(t){

	print(t)
	data_file <- list.files(path=input_dir, pattern=glue("_motif_enrichment_DNase_control_{t}.rds"))
	tissue_name <- sub("_motif_enrichment.*","",data_file)

	data_in <- lapply(data_file, function(x){
		readRDS(glue("{input_dir}/{x}"))
	})

	names(data_in) <- tissue_name

	data_OR <- lapply(data_in, function(x){
		tmp <- log2(x$odds_ratio)
		names(tmp) <- row.names(x)
		return(tmp)
	})

	return(data_OR)
})

#read in FDR
result_type_FDR <- lapply(region_type, function(t){

	print(t)
	data_file <- list.files(path=input_dir, pattern=glue("_motif_enrichment_DNase_control_{t}.rds"))
	tissue_name <- sub("_motif_enrichment.*","",data_file)

	data_in <- lapply(data_file, function(x){
		readRDS(glue("{input_dir}/{x}"))
	})

	names(data_in) <- tissue_name

	data_out <- lapply(data_in, function(x){
		tmp <- x$FDR
		names(tmp) <- row.names(x)
		return(tmp)
	})

	return(data_out)
})

##read in standard error of log odds ratio
result_type_SE <- lapply(region_type, function(t){

	print(t)
	data_file <- list.files(path=input_dir, pattern=glue("_motif_enrichment_DNase_control_{t}.rds"))
	tissue_name <- sub("_motif_enrichment.*","",data_file)

	data_in <- lapply(data_file, function(x){
		readRDS(glue("{input_dir}/{x}"))
	})

	names(data_in) <- tissue_name

	data_out <- lapply(data_in, function(x){
		tmp <- x$odds_ratio_se
		names(tmp) <- row.names(x)
		return(tmp)
	})

	return(data_out)
})

##fit linear regression model and use 75% prediction interval to identify differentially enriched motifs
level_cutoff <- 0.75
motif_residual <- lapply(names(result_type[[1]]), function(i){
	data_OR <- data.frame(dMML=result_type[[1]][[i]]/result_type_SE[[1]][[i]], dNME=result_type[[2]][[i]]/result_type_SE[[2]][[i]])
	data_FDR <- data.frame(dMML=result_type_FDR[[1]][[i]], dNME=result_type_FDR[[2]][[i]])

	model_fit <- lm(dNME ~ dMML, data_OR)
	pre_value <- predict(model_fit,interval="prediction",level=level_cutoff)

	residuals_out <- model_fit$residuals[(data_OR$dNME > pre_value[,3] & data_FDR$dNME < 0.1) |
																			(data_OR$dNME < pre_value[,2] & data_FDR$dMML < 0.1)]

	return(residuals_out)
})
names(motif_residual) <- names(result_type[[1]])

##compare with motifs that prefer high NME in human
motif_human_high_NME <- fread("../downstream/input/mouse_analysis/motif_prefer_high_NME_human.txt",header=FALSE) %>% as.data.frame()

##generate tables
sapply(names(motif_residual), function(x){
	motif_resudual_dNME <- motif_residual[[x]][motif_residual[[x]] > 0]
	motif_dNME <- sub(".*_","",names(motif_resudual_dNME))
	data_out <- data.frame(motif = names(motif_resudual_dNME),
												human_high_NME = motif_dNME %in% motif_human_high_NME[,1],
												residual = motif_resudual_dNME,
												log_OR_dNME = result_type[[2]][[x]][names(motif_resudual_dNME)],
												normalized_log_OR_dNME = result_type[[2]][[x]][names(motif_resudual_dNME)]/result_type_SE[[2]][[x]][names(motif_resudual_dNME)],
												FDR_dNME = result_type_FDR[[2]][[x]][names(motif_resudual_dNME)]
											)
	data_out <- data_out %>%
			arrange(desc(log_OR_dNME))

	write.csv(data_out,file=glue("{input_dir}/{x}_OR_residual_dNME.csv"), row.names=FALSE)

	motif_resudual_dMML <- motif_residual[[x]][motif_residual[[x]] < 0]
	motif_dMML <- sub(".*_","",names(motif_resudual_dMML))
	data_out <- data.frame(motif = names(motif_resudual_dMML),
												human_high_NME = motif_dMML %in% motif_human_high_NME[,1],
												residual = motif_resudual_dMML,
												log_OR_dMML = result_type[[1]][[x]][names(motif_resudual_dMML)],
												normalized_log_OR_dMML = result_type[[1]][[x]][names(motif_resudual_dMML)]/result_type_SE[[1]][[x]][names(motif_resudual_dMML)],
												FDR_dMML = result_type_FDR[[1]][[x]][names(motif_resudual_dMML)]
											)
	data_out <- data_out %>%
			arrange(desc(log_OR_dMML))

	write.csv(data_out,file=glue("{input_dir}/{x}_OR_residual_dMML.csv"), row.names=FALSE)
})

##motif enriched in MML only and NME only regions
motif_dMML <- lapply(motif_residual, function(i){
	names(i)[which(i < 0)]
})

motif_dMML_all <- Reduce(union, motif_dMML)

motif_dNME <- lapply(motif_residual, function(i){
		names(i)[which(i > 0)]
})

motif_dNME_all <- Reduce(union, motif_dNME)

motif_dNME_all <- sub(".*_","",motif_dNME_all)
dNME_hit <- length(which(motif_dNME_all %in% motif_human_high_NME[,1]))
dNME_off <- length(motif_dNME_all) - dNME_hit
overlap_motif_dNME <- motif_dNME_all[which(motif_dNME_all %in% motif_human_high_NME[,1])]
write.csv(overlap_motif_dNME, file="perfer_high_NME_overlap_motif_dNME.csv")

motif_dMML_all <- sub(".*_","",motif_dMML_all)
dMML_hit <- length(which(motif_dMML_all %in% motif_human_high_NME[,1]))
dMML_off <- length(motif_dMML_all) - dMML_hit
overlap_motif_dMML <- motif_dMML_all[which(motif_dMML_all %in% motif_human_high_NME[,1])]
write.csv(overlap_motif_dMML, file="perfer_high_NME_overlap_motif_dMML.csv")

##test for overlap between human and mouse results
test_table <- matrix(c(dNME_hit,dMML_hit,dNME_off,dMML_off),
											 nrow = 2,dimnames = list(c("dNME", "dMML"),c("hit", "off")))

test_result <- fisher.test(test_table, alternative = "greater")

# > test_table
#      hit off
# dNME  77 284
# dMML  56 295
# > test_result
#
# 	Fisher's Exact Test for Count Data
#
# data:  test_table
# p-value = 0.04043
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.01958     Inf
# sample estimates:
# odds ratio
#    1.42753
