rm(list=ls())
source("mainFunctions_sub.R")
theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()

# MAV vs NME -------------------------------------------------
NME_in_dt=readRDS(NME_mouse_MAV_fn)
NME_in_dt=NME_in_dt[(!is.na(hyper_var)&hyper_var!=-100)&!is.na(NME)]
NME_in_dt$Sample=NME_in_dt$stage
NME_in_dt$score=NME_in_dt$NME
dist_plot_run(NME_in_dt,theme_glob,ylab="NME",stat_in="hyper_var",dir=figure_path)
dist_plot_run(NME_in_dt,theme_glob,ylab="NME",stat_in="var",dir=figure_path)
dist_plot_run(NME_in_dt,theme_glob,ylab="NME",stat_in="mean",dir=figure_path)



#enhancer regions MAV vs NME
enhancer=readRDS(bin_enhancer_rds)
olap_enhancer=findOverlaps(convert_GR(NME_in_dt$region),enhancer)
NME_in_dt_enc=NME_in_dt[queryHits(olap_enhancer)]
NME_in_dt_enc$gene_enc=enhancer[subjectHits(olap_enhancer)]$`Target Gene`
NME_in_dt_enc$hyper_var_enc=-100
NME_in_dt_enc$var=-100
NME_in_dt_enc$mean=-100
for(st in unique(NME_in_dt_enc$stage)){
  tt1=proc.time()[[3]]
  if(file.exists(paste0(dir_scRNA_mouse,gsub('E|limb\\.|\\.all','',st),'.rds'))){
    scRNA_in=readRDS(paste0(dir_scRNA_mouse,gsub('E|limb\\.|\\.all','',st),'.rds'))
    scRNA_in=scRNA_in[rownames(scRNA_in)%in% unique(c(NME_in_dt[(stage==st)]$gene)),]
    if(nrow(scRNA_in)>0){
      #Add hypervar to TSS 
      NME_in_dt_enc[(stage==st)]$hyper_var_enc=scRNA_in[NME_in_dt_enc[(stage==st)]$gene_enc,"hypervar_logvar"]
      NME_in_dt_enc[(stage==st)]$var=scRNA_in[NME_in_dt_enc[(stage==st)]$gene,"var"]
      NME_in_dt_enc[(stage==st)]$mean=scRNA_in[NME_in_dt_enc[(stage==st)]$gene,"mean"]
      
    }
  }else{cat("File not exist for ",st,'\n')}
  cat('Finish processing ',sub('E','',st),'in: ',proc.time()[[3]]-tt1,'\n')
  
}
saveRDS(NME_in_dt_enc,NME_mouse_MAV_enhancer_fn)
NME_in_dt_enc_median_NME=NME_in_dt_enc[,list(NME=median(NME)),by=list(gene_enc,stage,hyper_var)]
NME_in_dt_enc_median_NME[gene_enc%in%names(table(NME_in_dt_enc_median_NME$gene_enc))[table(NME_in_dt_enc_median_NME$gene_enc)>=5],list(cor=cor(NME,hyper_var)),by=list(gene_enc)][order(cor,decreasing=T)]

NME_in_dt_enc[gene_enc%in%names(table(NME_in_dt_enc_median_NME$gene_enc))[table(NME_in_dt_enc_median_NME$gene_enc)>=5],list(cor=cor(NME,hyper_var)),by=list(gene_enc)][!is.na(cor)][order(cor,decreasing=T)]
cat(NME_in_dt_enc[gene_enc%in%names(table(NME_in_dt_enc_median_NME$gene_enc))[table(NME_in_dt_enc_median_NME$gene_enc)>=5],list(cor=cor(NME,hyper_var)),by=list(gene_enc)][cor>0.5]$gene_enc,sep='\n')

# motif preprocessing for Ken ----------------------------------------------------------
#See mouse_motif_processing.R


