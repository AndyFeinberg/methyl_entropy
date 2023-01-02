source('mainFunctions_sub.R')
tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
tissue_specific_GO=data.table(tissue_all=tissue_all,
                              GO_grepl=c("facial|odontogenesis|chondrocyte|myoblast|oteoblast",
                                         "axon|neuron|oligodendrocyte|synapse|ganglion|brain|forebrain|telencephalon|gyrus",
                                         "heart|cardiac|valve|aorta|ventricle|ventricular",
                                         "postsynapse|nervous|cerebellar|neural|hindbrain",
                                         "cartilage|bone|chondrocyte|epithelial",
                                         "hemopoiesis|leukocyte|hemopothesis|lymphocyte|monosaccharide",
                                         "glial|spinal cord|diencephalon|collateral sprouting|neuron")
)
tissue_out_filtered=readRDS(tissue_out_filtered_fn)
bin_enhancer=readRDS(bin_enhancer_rds)
tissue_out_filtered_enhancer=lapply(tissue_out_filtered,function(x) x[queryHits(findOverlaps(convert_GR(x$region),bin_enhancer))])
GO_out=readRDS(paste0(GO_01_dir,'GO_out_all_dMML_dNME_0rm_FC_N17_kmeans_10run_filtered_all_regions_01_enhancer.rds'))
UC_merge=readRDS(UC_merge_file)#Define all analyzed regions, were using UC_merge_max_loc_cluster01.rds,4626
uc_gr=lapply(UC_merge,function(x) rownames(x))
uc_gr=Reduce(intersect,uc_gr)
uc_gr=unique(uc_gr)
uc_enhancer_olap=findOverlaps(convert_GR(uc_gr),bin_enhancer)
uc_enhancer=uc_gr[queryHits(uc_enhancer_olap)]
percent_olap<-function(region_tissue,bin_enhancer,GO_tissue,tissue){
    #Annotate regions to genes 
    region_tissue_olap=findOverlaps(convert_GR(region_tissue$region),bin_enhancer)
    region_tissue_enhancer=region_tissue[queryHits(region_tissue_olap)]
    region_tissue_enhancer$gene="NA"
    region_tissue_enhancer$gene=bin_enhancer$`Target Gene`[subjectHits(region_tissue_olap)]
    GO_tissue=GO_out$all[[tissue]]
    GO_tissue_gene=do.call(rbind,lapply(GO_tissue,function(x) x$GO_out_cluster_all[FDR<=0.2&FC>1.5&grepl(tissue_specific_GO[tissue_all==tissue]$GO_grepl,Term)]))
    GO_gene=unique(unlist(strsplit(GO_tissue_gene$genes,';')))
    write.csv(GO_gene,paste0("../downstream/output/mouse_analysis/QC/ts_gene_",tissue,'.csv'))
    out=region_tissue_enhancer[,list(percent_GO=mean(gene%in%GO_gene)*100,total_GO=sum(gene%in%GO_gene),total=length(gene)),by=list(region_type)]
    out$expected_total_GO=out$total*out[region_type=="control"]$percent_GO/100
    out$tissue=tissue
    return(out)
}
tissue_UC_diff=lapply(tissue_all,function(tissue){
    #Find UC in the tissue
    region_tissue=tissue_out_filtered[[tissue]]
    #tissue_result=percent_olap(region_tissue,bin_enhancer,GO_tissue,tissue)
    #tissue_result$status="tissue"
    #Perfrom random control in all regions but not that tissue
    #region_non_tissue=do.call(rbind,tissue_out_filtered[names(tissue_out_filtered)!=tissue])
    #region_non_tissue=region_non_tissue[sample(uc_gr_tissue,nrow(region_tissue))]
    #region_non_tissue=uc_enhancer[sample(1:nrow(region_non_tissue),nrow(region_tissue))]
    uc_gr_tissue=which(!uc_enhancer%in%region_tissue$region)
    #region_non_tissue=uc_enhancer[sample(uc_gr_tissue,nrow(region_tissue))]
    region_non_tissue=data.table(region=uc_enhancer[uc_gr_tissue],region_type="control")
    region_all=rbind(region_tissue[,list(region,region_type)],region_non_tissue)
    non_tissue_result=percent_olap(region_all,bin_enhancer,GO_tissue,tissue)
    #non_tissue_result$status="control"
    return(non_tissue_result)
})
tissue_UC_diff=do.call(rbind,tissue_UC_diff)

tissue_UC_diff[,list(median_percent=median(percent_GO),mean_percent=mean(percent_GO),
                        sd_percent=sd(percent_GO),mean_region=mean(total_GO),mean_expected=mean(expected_total_GO)),by=list(region_type)]
write.csv(tissue_UC_diff,'../downstream/output/mouse_analysis/QC/ts_gene_tissue_number_summary.csv')
# theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=18),
#                                  axis.title.x=element_text(hjust=0.5,size=16,face="bold"),
#                                  axis.title.y=element_text(hjust=0.5,size=16,face="bold"),
#                                  axis.text.x=element_text(size=12),
#                                  axis.text.y=element_text(size=12))

# pdf("../downstream/output/mouse_analysis/QC/dMML_dNME_overlap_ts_gene_noeb.pdf",width=8,height=6)
# ggplot(tissue_UC_diff_bar,aes(x=region_type,fill=status,y=mean_percent))+
#             geom_bar(stat="identity", color="black", position=position_dodge()) +
#             #geom_errorbar(aes(ymin=median_percent-sd_percent, ymax=median_percent+sd_percent), width=.2, position=position_dodge(.9))+
#             xlab("region type")+ylab("percent overlap")+ggtitle("Percent region at the enhancers\nof genes with tissue-specific function")+
#             theme_glob+theme(legend.position="bottom")

# dev.off()
# pdf("../downstream/output/mouse_analysis/QC/dMML_dNME_overlap_ts_gene_eb.pdf",width=7,height=7)
# ggplot(tissue_UC_diff_bar,aes(x=region_type,fill=status,y=mean_percent))+
#             geom_bar(stat="identity", color="black", position=position_dodge()) +
#             geom_errorbar(aes(ymin=median_percent-sd_percent, ymax=median_percent+sd_percent), width=.2, position=position_dodge(.9))+
#             xlab("region type")+ylab("percent overlap")+ggtitle("Percent region at the enhancers\nof genes with tissue-specific function")+
#             theme_glob+theme(legend.position="bottom")

# dev.off()
#0.08
t.test(tissue_UC_diff[region_type=="Predominantly\nNME"]$percent_GO,tissue_UC_diff[region_type=="Predominantly\nMML"]$percent_GO,alternative="two.sided", paired = T)
#4.7*10^-3
t.test(tissue_UC_diff[region_type=="Predominantly\nNME"]$percent_GO,tissue_UC_diff[region_type=="control"]$percent_GO,alternative="two.sided", paired = T)
#2.6*10^-3
t.test(tissue_UC_diff[region_type=="Predominantly\nMML"]$percent_GO,tissue_UC_diff[region_type=="control"]$percent_GO,alternative="two.sided", paired = T)