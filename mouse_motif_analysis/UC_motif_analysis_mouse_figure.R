library(GenomicRanges)
library(dplyr)
library(glue)
library(ComplexHeatmap)
library(ggplot2)
library(patchwork)
library(data.table)
library(ggrepel)

##generate color
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color_2 <- gg_color_hue(2)

setwd("../downstream/output/mouse_analysis/UC_motif_enrichment")

region_type <- c("dMML", "dNME")
input_dir <- "all_region"

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

##read in FDR
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

##read in motifs that overlap with human analysis
select_motifs_dNME <- read.csv("perfer_high_NME_overlap_motif_dNME.csv")[,2] %>% as.character()
select_motifs_dMML <- read.csv("perfer_high_NME_overlap_motif_dMML.csv")[,2] %>% as.character()

select_motifs <- c(select_motifs_dNME, select_motifs_dMML)

motif_residual <- lapply(names(result_type[[1]]), function(i){

	data_OR <- data.frame(dMML=result_type[[1]][[i]]/result_type_SE[[1]][[i]], dNME=result_type[[2]][[i]]/result_type_SE[[2]][[i]])
	data_FDR <- data.frame(dMML=result_type_FDR[[1]][[i]], dNME=result_type_FDR[[2]][[i]])
	model_fit <- lm(dNME ~ dMML, data_OR)
	pre_value <- predict(model_fit,interval="prediction",level=0.75)

	residuals_out <- model_fit$residuals[(data_OR$dNME > pre_value[,3] & data_FDR$dNME < 0.05) |
																			(data_OR$dNME < pre_value[,2] & data_FDR$dMML < 0.05)]

	newx <- seq(min(data_OR$dMML),max(data_OR$dMML),by = 0.05)
	conf_interval <- predict(model_fit, newdata=data.frame(dMML=newx), interval="prediction",
													 level = 0.75)
	data_CI_upper <- data.frame(x=newx, y=conf_interval[,2])
	data_CI_lower <- data.frame(x=newx, y=conf_interval[,3])

	highlight_point_dNME <- data.frame(x=data_OR$dMML[(data_OR$dNME > pre_value[,3] & data_FDR$dNME < 0.05)],
																			y=data_OR$dNME[(data_OR$dNME > pre_value[,3] & data_FDR$dNME < 0.05)])

	highlight_point_dMML <- data.frame(x=data_OR$dMML[(data_OR$dNME < pre_value[,2] & data_FDR$dMML < 0.05)],
																			y=data_OR$dNME[(data_OR$dNME < pre_value[,2] & data_FDR$dMML < 0.05)])

	highlight_point_com <- rbind(highlight_point_dNME, highlight_point_dMML)

	motif_names <- c(row.names(data_OR)[(data_OR$dNME > pre_value[,3] & data_FDR$dNME < 0.05)],
									 row.names(data_OR)[(data_OR$dNME < pre_value[,2] & data_FDR$dMML < 0.05)])
	TF_names <- sub(".*_", "", motif_names)

	label_data_all <- data.frame(x = highlight_point_com$x,
														y = highlight_point_com$y,
														motif = motif_names,
														TF = TF_names,
														label = sub("\\(.*", "", TF_names),
														group = rep(c("dNME","dMML"),c(nrow(highlight_point_dNME),nrow(highlight_point_dMML))))

	label_data_mark <- label_data_all[label_data_all$TF %in% select_motifs,]

	##generate scatterplot
	p <- ggplot(data = data_OR, aes(x = dMML, y = dNME)) +
			geom_point() +
			xlab("Normalized Log2(odds-ratio) for dMML regions") +
			ylab("Normalized Log2(odds-ratio) for dNME regions") +
			ggtitle(i) +
			geom_abline(intercept = model_fit$coefficients[1], slope = model_fit$coefficients[2], col = "orange") +
			geom_line(data = data_CI_upper, aes(x = x, y = y), color = "grey", linetype = "dashed") +
			geom_line(data = data_CI_lower, aes(x = x, y = y), color = "grey", linetype = "dashed") +
			geom_point(data = label_data_all, aes(x = x, y = y, color = group)) +
			geom_label_repel(data = label_data_mark, aes(x = x, y = y, label = label, color = group), size = 5) +
			scale_colour_manual(values = c("dMML" = color_2[1],"dNME" = color_2[2])) +
			theme_bw() +
			theme(plot.title = element_text(size = 24, hjust = 0.5),
						axis.text = element_text(size = 15),
						axis.title = element_text(size = 15),
						legend.position = "none")
	ggsave(file=glue("{input_dir}/{i}_OR_compare_PI.svg"), plot=p, height = 5, width = 5)

	pdf(glue("{input_dir}/{i}_OR_compare_PI.pdf"), height = 5, width = 5)
	print(p)
	dev.off()

	return(residuals_out)
})
