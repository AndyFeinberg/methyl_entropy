source('mainFunctions_sub.R')
UC_in=readRDS(UC_in_matrix_cluster_file)
cut=0.1
#Non tissue-specific high UC
aid <- sapply(names(UC_in),function(i) {
 names(which(rowSums(UC_in[[i]] > cut) > 0))
  })  
aid_tb=table(table(unlist(aid)))
#Plotting histgram
aid_tb_dt=data.table(n_region=aid_tb)
colnames(aid_tb_dt)=c('n_tissue','n_region')
aid_tb_dt$n_region_pro=aid_tb_dt$n_region/sum(aid_tb_dt$n_region)
aid_tb_dt$n_tissue=factor(aid_tb_dt$n_tissue,levels=as.character(1:10))
pdf(paste0(figure_path,'number_tissue_pass_filter.pdf'))
ggplot(aid_tb_dt,aes(x=n_tissue,y=n_region_pro))+geom_bar(stat='identity', fill="steelblue")+
xlab('Number of tissue with regions having UC> 0.1')+ylab('Proportion of regions')+
geom_text(aes(label=round(n_region_pro,digits=2)),vjust=-1.6)
dev.off()
#Count overlap enhancers
aid_gr=lapply(aid,convert_GR,direction='GR')
enhancer_bin=readRDS(bin_enhancer_rds)
aid_gr_enhancer=lapply(aid_gr,function(x) subsetByOverlaps(x,enhancer_bin))
aid_gr_enhancer_dt=lapply(aid_gr_enhancer,function(x) convert_GR(x,direction="DT")$region)

aid_tb_enhancer=table(table(unlist(aid_gr_enhancer_dt)))
#  1     2     3     4     5     6     7
#26436 13077  8453  6645  7896  9589  7659
aid_tb_dt_enhancer=data.table(n_region=aid_tb_enhancer)
colnames(aid_tb_dt_enhancer)=c('n_tissue','n_region')
aid_tb_dt_enhancer$n_region_pro=aid_tb_dt_enhancer$n_region/sum(aid_tb_dt_enhancer$n_region)
aid_tb_dt_enhancer$n_tissue=factor(aid_tb_dt_enhancer$n_tissue,levels=as.character(1:10))
pdf(paste0(figure_path,'number_tissue_pass_filter_enhancer.pdf'))
ggplot(aid_tb_dt_enhancer,aes(x=n_tissue,y=n_region_pro))+geom_bar(stat='identity', fill="steelblue")+
xlab('Number of tissue with regions having UC> 0.1')+ylab('Proportion of regions')+
geom_text(aes(label=round(n_region_pro,digits=2)),vjust=-1.6)
dev.off()

#Generate heatmap using non-ts cluster
#Generate cluster
for (seed in 1:10) {
  cut <-0.1
  aid <- sapply(names(UC_in),function(i) {
    names(which(rowSums(UC_in[[i]] > cut) > 0))
  })  
   timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
  d2 <- sapply(names(UC_in),function(i) {
    ###### this is the line to ensure it's > cut in only one tissue
    #sid <- setdiff(aid[[i]],unlist(aid[names(aid)!=i]))
    sid=aid[[i]]
    i <- UC_in[[i]]
    i <- i[sid,]
    i <- i[,colnames(i) %in% timeorder]
    scalematrix <- function(data) {
      cm <- rowMeans(data)
      csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
      (data - cm) / csd
    }
    i <- scalematrix(i)
    i <- i[complete.cases(i),]
    set.seed(seed)
    clu <- kmeans(i,10,iter.max = 10000)$cluster
    n <- names(clu)
    clum <- rowsum(i,clu)/as.vector(table(clu))
    maxp <- apply(clum,1,function(i) {
      i[i < 0] <- 0
      i <- i-min(i)
      i <- i/max(i)
      sum(i*c(1:length(i)))/sum(i)
    })
    clu <- rank(maxp)[clu]
    names(clu) <- n
    clu
  })
  saveRDS(d2,file=paste0('../downstream/output/mouse_analysis/non_ts_clustering/uc_',cut,'_',seed,'.rds'))
}
UC_in=readRDS(UC_merge_max_loc_file)
cutoff_char='01'
cluster_assignment('../downstream/output/mouse_analysis/non_ts_clustering/',
                    '../downstream/output/mouse_analysis/non_ts_clustering/cluster_assigned/',
                    UC_merge=UC_in,
                    cutoffs=0.1,
                    cluster_region_out_fn=paste0('../downstream/output/mouse_analysis/non_ts_clustering/cluster_assginment_filtered_',cutoff_char,'.rds'),
                    figure_name='../downstream/output/mouse_analysis/non_ts_clustering/heatmap_non_ts.png',
                    figure_width=2000,figure_height=20000,res=200)
# clu=readRDS( paste0('../downstream/output/mouse_analysis/non_ts_clustering/cluster_assginment_filtered_',cutoff_char,'.rds'))
# plot_heatmap_cluster(UC_in,clu,
#                     figure_name='../downstream/output/mouse_analysis/non_ts_clustering/heatmap_non_ts.png',
#                     figure_width=2000,figure_height=20000,res=200)