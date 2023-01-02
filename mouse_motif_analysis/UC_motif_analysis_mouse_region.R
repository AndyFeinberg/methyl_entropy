library(data.table)
library(dplyr)
library(glue)
library(GenomicRanges)
options(scipen=999)

setwd("../downstream/input/mouse_analysis")

##read in all regions and retain enhancer regions
data_list <- readRDS("../downstream/input/mouse_analysis/tissue_out_N17_kmeans_10run_filtered_all_region.rds")
tissue_name <- names(data_list)

cluster_data <- lapply(data_list, function(x){
	x <- x[x$enhancer == TRUE, c("region","cluster","region_type")]
})
names(cluster_data) <- tissue_name

cluster_gr <- lapply(cluster_data, function(data_in){
	peak_data <- data_in$region %>% sapply(.,function(x)strsplit(x,paste("-",":",sep="|"))[[1]]) %>% t()
	peak_gr <- GRanges(seqnames=peak_data[,1],ranges=IRanges(start=as.numeric(peak_data[,2]),end=as.numeric(peak_data[,3])),cluster=data_in$cluster, region_type=data_in$region_type)
})

saveRDS(cluster_gr, file="tissue_region_enhancer_gr.rds")

##get control regions using cisgenome for each tissue
cluster_gr_control <- lapply(1:length(cluster_gr), function(i){
	gr_cod <- data.frame(cluster_gr[[i]])[,c(1:3)]
	gr_cod$strand <- "+"
	write.table(gr_cod, file=glue("{names(cluster_gr)[i]}.cod"), col.names=FALSE, sep="\t", quote=FALSE)
	system(glue("~/cisgenome_project/bin/refgene_getmatchedcontrol -d /users/wzhou14/cisgenome_project/mm10/annotation/refFlat_sorted.txt -dt 1 -s mouse -c /users/wzhou14/cisgenome_project/mm10/chrlen.txt -i {names(cluster_gr)[i]}.cod -o {names(cluster_gr)[i]}_control.cod -n 1 -l 250 -nr 1"))
	gr_control <- fread(glue("{names(cluster_gr)[i]}_control.cod"))
	gr_control <- makeGRangesFromDataFrame(gr_control[,-1])
	return(gr_control)
})

names(cluster_gr_control) <- names(cluster_gr)

saveRDS(cluster_gr_control, file="tissue_region_control_enhancer_gr.rds")
