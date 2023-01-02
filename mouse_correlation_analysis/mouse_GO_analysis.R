rm(list=ls())
source('mainFunctions_sub.R')
# GO analysis -------------------------------------------------------------
#Define input parameters
UC_merge=readRDS(UC_merge_file)#Define all analyzed regions, were using UC_merge_max_loc_cluster01.rds,4626
cutoff_fn='01'
#Runnning
tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
#prepare enhancer background gene list
uc_gr=lapply(UC_merge,function(x) rownames(x))
uc_gr=Reduce(intersect,uc_gr)
uc_gr=convert_GR(uc_gr)
enhancer=readRDS(bin_enhancer_rds)#21441
enhancer_bg=subsetByOverlaps(enhancer,uc_gr)
bg_enhancer=unique(enhancer_bg$`Target Gene`)
#Prepare promoter background gene
tss=get_mm10_tss()

bg_promoter=names(subsetByOverlaps(tss,uc_gr,maxgap = 2000))


# GO run for promoters, enhancers and different catogries -----------------
for(enc_type in c("enhancer","promoter")){
  enc_type="enhancer"
  GO_out_all=list()
  for(region_type in c("all","NME only","Neither","Both","MML only")){
    GO_out_all[[region_type]]=list()
    for(ts in tissue_all){
      if(enc_type=="enhancer"){
        bg=bg_enhancer
      }else 
        if(enc_type=="promoter"){
          bg=bg_promoter
          
        }
     
      GO_out_all[[region_type]][[ts]]=GO_run_tissue(ts,dir_out_cluster01,enc_type=enc_type,region_type_sel=region_type,bg=bg,DNase=F)
      GO_out_all[[region_type]][[ts]]=lapply(GO_out_all[[region_type]][[ts]],function(x){
        return(list(GO_out_cluster_all=x$GO_out_cluster_all,
                    csv_in_ts_clu=cbind(x$csv_in_ts_clu,as.data.table(UC_merge[[ts]][x$csv_in_ts_clu$region,!grepl('max',colnames(UC_merge[[ts]]))]))))
       
      })
      
    }
  }
  
  saveRDS(GO_out_all,paste0(GO_01_dir,'GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_',cutoff_fn,'_',enc_type,'.rds'))
}

# Plot heatmaps -----------------------------------------------------------


GO_out_all=readRDS(GO_01_enhancer_fn)
tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
#Plot all terms in a single plot
for(region_type in names(GO_out_all)){
  plot_GO_heatmap_all(tissue_all,GO_out_all[[region_type]],region_type=region_type,enc_type="enhancer",ptcount=0,FDR_cutoff=0.2,
                      dir_plot=GO_01_dir)
  
}
enc_type='promoter'
cutoff_fn='01'
GO_out_all=readRDS(GO_01_promoter_fn)
for(region_type in names(GO_out_all)){
  plot_GO_heatmap_all(tissue_all,GO_out_all[[region_type]],region_type=region_type,enc_type="promoter",ptcount=0,FDR_cutoff=0.2,
                      dir_plot=GO_01_dir)
  
}


# Write csv output --------------------------------------------------------

chrs <- names(Mmusculus)[1:21]#2276796
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Mmusculus[[x]])))
cpgr <- do.call(c, lapply(1:21, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs[[x]], width = 1)))) #use first location
GO_out_all=readRDS(paste0(GO_01_dir,'GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_',cutoff_fn,'_',enc_type,'.rds'))
GO_sheets(GO_out_all,"enhancer",dMML_cor=dMML_cor,dNME_cor=dNME_cor,mm10_CpG=cpgr,FDR_cutoff = 0.2,out_dir=paste0(GO_01_dir,'GO_sheets/'))
enc_type="promoter"
GO_out_all=readRDS(paste0(GO_01_dir,'GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_',cutoff_fn,'_',enc_type,'.rds'))
GO_sheets(GO_out_all,"promoter",dMML_cor=dMML_cor,dNME_cor=dNME_cor,mm10_CpG=cpgr,FDR_cutoff = 0.2,out_dir=paste0(GO_01_dir,'GO_sheets/'))


