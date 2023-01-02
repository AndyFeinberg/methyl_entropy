source('mainFunctions_sub.R')


plot_out_cor=list()
for(cutoff_char in  c("0025","005","015","02")){
    output_dir=paste0('../downstream/output/mouse_analysis/cutoff_testing/',cutoff_char,'/')
    plot_out_all=paste0(output_dir,'tissue_out_N17_kmeans_10run_filtered_all_region_',cutoff_char,'.rds')
   plot_cutoff=plot_correlation(readRDS(plot_out_all),pdf_fn=paste0(output_dir,cutoff_char,'_correlation_main.pdf'),plot_pdf=F)+
                ggtitle(paste0("cutoff:",cutoff_char))+theme(plot.title = element_text(size=36,hjust=0.5))
    plot_out_cor[[cutoff_char]]=plot_cutoff

                
}
ggsave(
    filename='../downstream/output/mouse_analysis/cutoff_testing/correlation_fig2.pdf',
    plot=ggarrange(plotlist=plot_out_cor, ncol=2, nrow=2, common.legend = TRUE, legend="bottom"),
   width = 21,
   height =28,
)

#heatmap for values
d=readRDS(UC_merge_max_loc_file)
for(cutoff_char in  c("0025","015")){
    output_dir=paste0('../downstream/output/mouse_analysis/cutoff_testing/',cutoff_char,'/')
    cluster_assigned_dir=paste0(output_dir,'cluster_assigned/')
    figure_name=paste0(figure_path,'all_sc_N17_ft_kmeans_10run_filtered_all',cutoff_char,'.png')
    clu=readRDS(paste0(cluster_assigned_dir,'cluster_assginment_filtered_',cutoff_char,'.rds'))
    #0.2, 0.3, 0.005 using width=1800, height=2000
    #0.15: width=1800, height=2500, res=300
    #0.025: width=1800, height=4000, res=300
    if(cutoff_char =="0025"){
        figure_height=4000
        res=300
    }else if (cutoff_char=="015"){
        figure_height=2500
        res=300
    }else if(cutoff_char %in% c("0.005","0.2")){
        figure_height=2000
        res=200
    }
    plot_heatmap_cluster(d,clu,figure_name,figure_width=1800,figure_height=figure_height,res=res)
}
#Get quantile of each cutoffs
UC=readRDS(UC_merge_file)
UC=lapply(UC,function(x) x[,grepl("UC-",colnames(x))])
quantile(fastDoCall('c',lapply(UC,as.vector)),c(0.7,0.85,0.95,0.97,0.99))
#        70%        85%        95%        97%        99% 
# 0.02245275 0.04461095 0.12623131 0.16390048 0.2241757
UC_ecdf=ecdf(fastDoCall('c',lapply(UC,as.vector)))
saveRDS(UC_ecdf,'../downstream/output/mouse_analysis/cutoff_testing/UC_ecdf.rds')
UC_ecdf_rowMax=ecdf(fastDoCall('c',lapply(UC,rowMaxs)))
saveRDS(UC_ecdf_rowMax,'../downstream/output/mouse_analysis/cutoff_testing/UC_ecdf_rowMax.rds')

 
#Plotting the UC01 ones
figure_name=paste0(figure_path,'all_sc_N17_ft_kmeans_10run_filtered_all',gsub('.','',cutoffs),'.tiff')
clu=readRDS(cluster_01_region_out_fn)
plot_heatmap_cluster(d,clu,figure_name,figure_width=1800,figure_height=2000,res=200)
