library(GenomicRanges)
library(dplyr)
library(glue)
library(ComplexHeatmap)
library(ggplot2)
library(patchwork)
library(ggsignif)
library(matrixStats)
source('mainFunctions_sub.R')


##read in NME data
data_median_com_DNase <- readRDS(hg19_DNase_NME_fn)
data_median_com_control <- readRDS(hg19_control_NME_fn)

##get sample order based on column mean
data_median_com_DNase_colmean <- colMeans(data_median_com_DNase, na.rm=TRUE)
col_order <- sort(data_median_com_DNase_colmean, decreasing = FALSE, index.return=TRUE)$x %>% names()

data_median_com_filter <- cbind(data_median_com_control[,col_order], data_median_com_DNase[,col_order])

##get motif order based on the row mean difference between regulatory and non-regulatory DNA regions
motif_NME_mean_DNase_combine_ave <- rowMeans(data_median_com_DNase, na.rm=TRUE)
motif_NME_mean_control_combine_ave <- rowMeans(data_median_com_control, na.rm=TRUE)

motif_order <- sort(motif_NME_mean_DNase_combine_ave - motif_NME_mean_control_combine_ave, decreasing = TRUE, index.return=TRUE)$x %>% names()

##standardize the NME values across regulatory and non-regulatory DNA regions for each motif
data_median_com_filter_sd <- apply(data_median_com_filter, 1, scale) %>% t()
dimnames(data_median_com_filter_sd) <- dimnames(data_median_com_filter)
data_median_com_filter_sd <- data_median_com_filter_sd[motif_order,]

##generate heatmap for all motifs
region_group <- rep(c("Non-regulatory DNA","Regulatory DNA"),each=49)

row_labels <- row.names(data_median_com_filter_sd) %>% sub(".*_","",.)
names(row_labels) <- row.names(data_median_com_filter_sd)

ht <- Heatmap(data_median_com_filter_sd, name = "NME",
							cluster_rows = FALSE,
							cluster_columns = FALSE,
							column_split = factor(region_group,levels=c("Non-regulatory DNA","Regulatory DNA")),
							column_gap = unit(5, "mm"),
              show_row_names = TRUE,
							row_labels = row_labels,
							show_column_names = TRUE,
							column_names_rot = 90,
              row_title_gp = gpar(fontsize = 40), row_names_gp = gpar(fontsize = 6),
							column_title_gp = gpar(fontsize = 40), column_names_gp = gpar(fontsize = 8),
							heatmap_legend_param = list(title_gp = gpar(fontsize = 30), labels_gp = gpar(fontsize = 24),
							legend_height = unit(3, 'in'), grid_width = unit(0.5, 'in'))

)

pdf(paste0(figure_path,"motif_NME_heatmap_median_human_DNase_control_agnostic.pdf"), width=20, height=40, family='ArialMT', useDingbats=FALSE)
draw(ht)
dev.off()

##generate heatmap for the top 50 motifs
data_median_com_filter_top50 <- cbind(data_median_com_control[motif_order[1:50],col_order], data_median_com_DNase[motif_order[1:50],col_order])

data_median_com_filter_top50_sd <- apply(data_median_com_filter_top50, 1, scale) %>% t()
dimnames(data_median_com_filter_top50_sd) <- dimnames(data_median_com_filter_top50)

row_labels <- row.names(data_median_com_filter_top50_sd) %>% sub(".*_","",.)
names(row_labels) <- row.names(data_median_com_filter_top50_sd)

ht <- Heatmap(data_median_com_filter_top50_sd, name = "NME",
							cluster_rows = FALSE,
							cluster_columns = FALSE,
							column_split = factor(region_group,levels=c("Non-regulatory DNA","Regulatory DNA")),
							column_gap = unit(5, "mm"),
              show_row_names=TRUE,
							row_labels = row_labels,
							show_column_names=TRUE,
							column_names_rot = 90,
              row_title_gp = gpar(fontsize = 30), row_names_gp = gpar(fontsize = 14),
							column_title_gp = gpar(fontsize = 30), column_names_gp = gpar(fontsize = 12)
)

pdf(paste0(figure_path,"motif_NME_heatmap_median_human_DNase_control_agnostic_top50.pdf"), width=16, height=15, family='ArialMT', useDingbats=FALSE)
draw(ht,padding = unit(c(40, 2, 2, 2), "mm"))
dev.off()

##perform differential test between motif sites at non-regulatory DNA and regulatory DNA regions for each motif
data_all <- data.frame(group = rep(c("Non-regulatory DNA","Regulatory DNA"), each = 49), t(data_median_com_filter))

data_test_p <- lapply(colnames(data_all[,-1]), function(t){
  data_sub <- data_all[,c("group",t)]
	colnames(data_sub) <- c("group","NME")
	p_value <- wilcox.test(data_sub[data_sub$group=="Non-regulatory DNA","NME"],data_sub[data_sub$group=="Regulatory DNA","NME"],paired = TRUE)$p.value
})

data_test_p <- Reduce(c, data_test_p)
data_test_FDR <- p.adjust(data_test_p, method="BH")

data_all_stat <- colMeans(data_all[data_all$group == "Regulatory DNA", -1], na.rm=T) - colMeans(data_all[data_all$group == "Non-regulatory DNA", -1], na.rm=T)
data_all_stat <- data.frame(motif = names(data_all_stat),TF=sapply(names(data_all_stat), function(x) strsplit(x,"_")[[1]][2]),
															diff=data_all_stat, p_value=data_test_p, FDR=data_test_FDR)

data_all_stat <- data_all_stat %>%
	mutate(mark = case_when(
		FDR < 0.001 ~ "***",
		FDR < 0.01 ~ "**",
		FDR < 0.05 ~ "*",
		TRUE ~ "NS"
	))

write.csv(data_all_stat, file=paste0(figure_path,"NME_diff_test_TFBS_human.csv"))

##generate violin plot for the top top 50 motifs
data_plot_norm <- data_median_com_filter[motif_order,]
data_all_plot <- data.frame(group = factor(rep(c("Non-regulatory","Regulatory"), each = 49),
levels=c("Non-regulatory","Regulatory")), t(data_plot_norm[1:50,]))

plot_list <- lapply(colnames(data_all_plot[,-1]), function(t){

  data_sub <- data_all_plot[,c("group",t)]
	colnames(data_sub) <- c("group","NME")

	p <- ggplot(data_sub, aes(x=group, y=NME, fill=group)) +
		geom_violin(scale = 'width') +
		geom_point(pch = 21, position = position_jitterdodge()) +
		geom_signif(annotations = data_all_stat$mark[data_all_stat$motif == t], y_position = 1, xmin="Non-regulatory", xmax="Regulatory") +
		ylim(c(0, 1.15)) +
		ggtitle(t) +
		theme_bw() +
		theme(plot.title = element_text(size = 14, hjust = 0.5),
					axis.text = element_text(size = 12),
					axis.title = element_blank())
})

pdf(paste0(figure_path,"motif_NME_violin_median_human_DNase_control_agnostic_top50.pdf"), width=20, height=38, family='ArialMT', useDingbats=FALSE)
wrap_plots(plot_list, ncol=5, guides="collect")
dev.off()
