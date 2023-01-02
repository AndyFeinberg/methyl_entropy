# nolint start
args = commandArgs(trailingOnly=TRUE)
cutoff=args[1]
cat("Processing:",cutoff,'\n')
source('mainFunctions_sub.R')
#We have already done the clustering for each cutoff, we need to run cluster assignment and GO and correlation analysis# nolint
#Input folder 
cutoff_char=as.character(gsub('\\.','',cutoff))
cluster_dir_in=paste0('../downstream/input/mouse_analysis/clustering/tissue_specific/uc_',cutoff_char,'/')
#creating output folder
output_dir=paste0('../downstream/output/mouse_analysis/cutoff_testing/',cutoff_char,'/')
if(!dir.exists(output_dir)){dir.create(output_dir)}
cluster_assigned_dir=paste0(output_dir,'cluster_assigned/')
if(!dir.exists(cluster_assigned_dir)){dir.create(cluster_assigned_dir)}
  UC_merge=readRDS(UC_merge_max_loc_file)
#Clustering
cat("Start cluster generation\n")
tt1=proc.time()[[3]]
figure_name=paste0(figure_path,'all_sc_N17_ft_kmeans_10run_filtered_all',cutoff_char,'.png')
cluster_assignment(cluster_dir_in,cluster_assigned_dir,cutoffs=cutoff,
            cluster_region_out_fn=paste0(cluster_assigned_dir,'cluster_assginment_filtered_',cutoff_char,'.rds'),figure_name=figure_name,UC_merge)

#Reading in all correlation
dir_in_cor_Jason='../downstream/input/mouse_analysis/correlation_analysis/all_regions/'
dMML_cor=readRDS(dmml_cor_file)
dNME_cor=readRDS(dnme_cor_file)
dmml_perm=readRDS(paste0(dir_in_cor_Jason,'fullpermudmmlcor_YQ.rds'))
dnme_perm=readRDS(paste0(dir_in_cor_Jason,'fullpermudnmecor_YQ.rds'))
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"), # nolint
                                 axis.text.x=element_text(size=16),
                                 axis.text.y=element_text(size=16))
cat("Finish cluster generation in:", proc.time()[[3]]-tt1,"\n")
tt1=proc.time()[[3]]
#Filtered result
cat("Start correlation analysis\n")
cor_dt_pre_process_fn=paste0(output_dir,'correlation_dt_N17_kmeans_10run_filtered_all_regions_',cutoff_char,'.rds')
if(!file.exists(cor_dt_pre_process_fn)){
  cor_dt_filtered=lapply(names(dNME_cor),cor_dt_preprocessing,dMML_cor=dMML_cor,
                         dNME_cor=dNME_cor,dmml_perm=dmml_perm,dnme_perm=dnme_perm,
                         filtered=TRUE,folder_input=cluster_assigned_dir)
  names(cor_dt_filtered)=names(dNME_cor)
  saveRDS(cor_dt_filtered,cor_dt_pre_process_fn)
}
cor_dt_filtered=readRDS(cor_dt_pre_process_fn)
#Plot the density for each one
library(ggpubr)
fig_out_dir=paste0(output_dir,'correlation_figure/')
cor_dt_pre_process_fn=paste0(output_dir,'tissue_out_N17_kmeans_10run_filtered_all_region_',cutoff_char,'.rds')
if(!dir.exists(fig_out_dir)){dir.create(fig_out_dir)}
tissue_out_filtered=lapply(names(cor_dt_filtered),correlation_processing,cor_dt=cor_dt_filtered,filtered=T,
                           dir_figure=fig_out_dir)
names(tissue_out_filtered)=names(cor_dt_filtered)

saveRDS(tissue_out_filtered,cor_dt_pre_process_fn)
plot_correlation(tissue_out_filtered,pdf_fn=paste0(output_dir,cutoff_char,'_correlation_main.pdf'))
cat("Finish correlation analysis in:", proc.time()[[3]]-tt1,"\n")
tt1=proc.time()[[3]]
cat("Starting GO analysis\n")
#GO analysis for enhancers
UC_merge=readRDS(UC_merge_file)#Define all analyzed regions, were using UC_merge_max_loc_cluster01.rds,4626
cutoff_fn=cutoff_char
#Runnning
tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
#prepare enhancer background gene list
uc_gr=lapply(UC_merge,function(x) rownames(x))
uc_gr=Reduce(intersect,uc_gr)
uc_gr=convert_GR(uc_gr)
enhancer=readRDS(bin_enhancer_rds)#21441
enhancer_bg=subsetByOverlaps(enhancer,uc_gr)
bg_enhancer=unique(enhancer_bg$`Target Gene`)

# GO run for promoters, enhancers and different catogries -----------------
enc_type="enhancer"
GO_out_all=list()
region_type ="all"
for(ts in tissue_all){
    GO_out_all[[ts]]=GO_run_tissue(ts,cluster_assigned_dir,enc_type=enc_type,region_type_sel=region_type,bg=bg_enhancer)
    GO_out_all[[ts]]=lapply(GO_out_all[[ts]],function(x){
                                                                            return(list(GO_out_cluster_all=x$GO_out_cluster_all,
                                                                            csv_in_ts_clu=cbind(x$csv_in_ts_clu,
                                                                            as.data.table(UC_merge[[ts]][x$csv_in_ts_clu$region,
                                                                            !grepl('max',colnames(UC_merge[[ts]]))]))))
        
        })
}

saveRDS(GO_out_all,paste0(output_dir,'GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_',cutoff_fn,'_enhancer.rds'))
GO_out_all=readRDS(paste0(output_dir,'GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_',cutoff_fn,'_enhancer.rds'))
tissue_select=names(which(table(select_top_GO(GO_out_all,names(GO_out_all),ptcount=0,FDR_cutoff=0.2,FC_cutoff=1.5)[[1]][FDR<=0.2&FC>=1.5]$tissue)>0))
if(!is.null(tissue_select)){
plot_GO_heatmap_all(tissue_select,GO_out_all,region_type="all",enc_type="enhancer",ptcount=0,FDR_cutoff=0.2,
                      dir_plot=output_dir)
}else{cat("No significant GO for cutoff:",cutoff,'\n')}
cat("Finish GO analysis in:", proc.time()[[3]]-tt1,"\n")
tt1=proc.time()[[3]]
# nolint end