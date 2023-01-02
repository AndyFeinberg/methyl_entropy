library(GenomicRanges)
library(dplyr)
library(glue)
library(ComplexHeatmap)
library(ggplot2)
library(patchwork)
library(ggsignif)
library(matrixStats)

##function for color generation
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

setwd("../downstream/output/mouse_analysis/TFBS_analysis")

##read in NME data
data_median_com_DNase <- readRDS("mouse_motif_DNase_median.rds")
data_median_com_control <- readRDS("mouse_motif_control_median.rds")

##exclude P0 time point and tissues not used
data_median_com_DNase <- data_median_com_DNase[,!grepl("P0|Lung|NT|intestine|kidney|stomach",colnames(data_median_com_DNase))]
data_median_com_control <- data_median_com_control[,!grepl("P0|Lung|NT|intestine|kidney|stomach",colnames(data_median_com_control))]

data_median_com_filter <- cbind(data_median_com_control,
																data_median_com_DNase)

##get motif order based on the row mean difference between regulatory and non-regulatory DNA regions
motif_NME_mean_DNase_combine_ave <- rowMeans(data_median_com_DNase, na.rm=TRUE)
motif_NME_mean_control_combine_ave <- rowMeans(data_median_com_control, na.rm=TRUE)

motif_order <- sort(motif_NME_mean_DNase_combine_ave - motif_NME_mean_control_combine_ave, decreasing = TRUE, index.return=TRUE)$x %>% names()

##standardize the NME values across regulatory and non-regulatory DNA regions for each motif
data_median_com_filter_sd <- apply(data_median_com_filter, 1, scale) %>% t()
dimnames(data_median_com_filter_sd) <- dimnames(data_median_com_filter)
colnames(data_median_com_filter_sd) <- sub("-all","",colnames(data_median_com_filter_sd))
data_median_com_filter_sd <- data_median_com_filter_sd[motif_order,]

##generate heatmap for all motifs
tissue_group <- colnames(data_median_com_filter) %>% sub("-.*","",.)
region_group <- rep(c("non-regulatory","regulatory"),each=46)

cols_tissue = gg_color_hue(length(unique(tissue_group)))
names(cols_tissue) <- unique(tissue_group)

column_bar1 <-  HeatmapAnnotation(tissue = tissue_group[region_group == "non-regulatory"],
col = list(tissue = cols_tissue),
show_annotation_name = FALSE,
annotation_legend_param = list(tissue = list(title_gp = gpar(fontsize = 30), labels_gp = gpar(fontsize = 20),
grid_height = unit(0.5, 'in'), grid_width = unit(0.5, 'in')))
)

ht1 <- Heatmap(data_median_com_filter_sd[,region_group == "non-regulatory"], name = "NME", top_annotation = column_bar1,
							cluster_rows = FALSE,
							cluster_columns = FALSE,
							column_split = factor(tissue_group[region_group == "non-regulatory"]),
							column_gap = unit(2, "mm"),
              show_row_names=TRUE,
							show_column_names=TRUE,
							column_names_rot = 90,
							column_title = "non-regulatory",
              row_title_gp = gpar(fontsize = 40), row_names_gp = gpar(fontsize = 6),
							column_title_gp = gpar(fontsize = 40), column_names_gp = gpar(fontsize = 8),
							heatmap_legend_param = list(title_gp = gpar(fontsize = 30), labels_gp = gpar(fontsize = 24),
							legend_height = unit(3, 'in'), grid_width = unit(0.5, 'in'))
)

column_bar2 <-  HeatmapAnnotation(tissue = tissue_group[region_group == "regulatory"],
col = list(tissue = cols_tissue),
annotation_legend_param = list(tissue = list(title_gp = gpar(fontsize = 30), labels_gp = gpar(fontsize = 24),
grid_height = unit(0.5, 'in'), grid_width = unit(0.5, 'in')))
)

ht2 <- Heatmap(data_median_com_filter_sd[, region_group == "regulatory"], name = "NME", top_annotation = column_bar2,
							cluster_rows = FALSE,
							cluster_columns = FALSE,
							column_split = factor(tissue_group[region_group == "regulatory"]),
							column_gap = unit(2, "mm"),
              show_row_names=TRUE,
							show_column_names=TRUE,
							column_names_rot = 90,
							column_title = "regulatory",
              row_title_gp = gpar(fontsize = 40), row_names_gp = gpar(fontsize = 6),
							column_title_gp = gpar(fontsize = 40), column_names_gp = gpar(fontsize = 8),
							heatmap_legend_param = list(title_gp = gpar(fontsize = 30), labels_gp = gpar(fontsize = 24),
							legend_height = unit(3, 'in'), grid_width = unit(0.5, 'in'))
)

pdf("motif_NME_heatmap_median_mouse_DNase_control.pdf", width=24, height=44, family='ArialMT', useDingbats=FALSE)
draw(ht1+ht2, ht_gap = unit(1, "cm"), column_title_gp = gpar(fontsize = 40))
dev.off()

##perform differential test between motif sites at non-regulatory DNA and regulatory DNA regions for each motif
data_all <- data.frame(group = rep(c("non-regulatory","regulatory"), each = 46), t(data_median_com_filter))

data_test_p <- lapply(colnames(data_all[,-1]), function(t){
  data_sub <- data_all[,c("group",t)]
	colnames(data_sub) <- c("group","NME")
	p_value <- wilcox.test(data_sub[data_sub$group=="non-regulatory","NME"],data_sub[data_sub$group=="regulatory","NME"],paired = TRUE)$p.value
})

data_test_p <- Reduce(c, data_test_p)
data_test_FDR <- p.adjust(data_test_p, method="BH")

data_all_stat <- colMeans(data_all[data_all$group == "regulatory", -1], na.rm=T) - colMeans(data_all[data_all$group == "non-regulatory", -1], na.rm=T)
data_all_stat <- data.frame(motif = names(data_all_stat),TF=sapply(names(data_all_stat), function(x) strsplit(x,"_")[[1]][2]),
															diff=data_all_stat, p_value=data_test_p, FDR=data_test_FDR)

data_all_stat <- data_all_stat %>%
	mutate(mark = case_when(
		FDR < 0.001 ~ "***",
		FDR < 0.01 ~ "**",
		FDR < 0.05 ~ "*",
		TRUE ~ "NS"
	))

write.csv(data_all_stat, file="diff_test_TFBS_mouse.csv")
