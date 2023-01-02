rm(list=ls())
source("mainFunctions_sub.R")
#Define ggplot theme

theme_glob=theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme_classic()
# add hypervar to allele-agnostic data ------------------------------------


NME_in=readRDS(NME_agnostic_file)
MML_in=readRDS(MML_agnostic_file)
genomic_features=readRDS(genomic_features_file)
GR_calc=data.frame()
scRNA_result=data.frame()
MML_hypervar_calc=GRanges()
NME_hypervar_calc=GRanges()
#For the mean expression: log2 or not?
for (sp in unique(NME_in$Sample)){
  
  hyper_var_file=unlist(strsplit(unique(NME_in$hyper_var_fn[NME_in$Sample==sp]),';'))
  hyper_var_file=gsub(scRNA_dir,scRNA_dir,hyper_var_file)
  cat('Processing',sp,'\n')
  if(all(file.exists(hyper_var_file))){
    
    sp_hyper_var=read_hypervar(hyper_var_file)
    print(head(sp_hyper_var))
    NME_hypervar_calc=c(NME_hypervar_calc,dist_plot_calc(NME_in[NME_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    MML_hypervar_calc=c(MML_hypervar_calc,dist_plot_calc(MML_in[MML_in$Sample==sp],sp_hyper_var,
                                                         genomic_features))
    
    
  }else{cat("file not exist for:",sp,'\n')}
}

#NME: 53476752 check
#MML: check: 53478342
saveRDS(list(NME_hypervar_calc=NME_hypervar_calc,
             MML_hypervar_calc=MML_hypervar_calc),
             paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV.rds'))
# Find number of overlapped regions ---------------------------------------

hyper_var_all=readRDS(paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV.rds'))#cor=0.211722 
hyper_var_all=lapply(hyper_var_all,function(x) x[x$N>=2])
hyper_var_all=lapply(hyper_var_all,convert_GR,direction='DT')
#Figure 3C and D in different context
#0.2324964
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="hypervar_logvar",dir=NME_MAV_human_out_dir)
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="mean",dir=NME_MAV_human_out_dir)
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="var",dir=NME_MAV_human_out_dir)
hyper_var_all$NME_hypervar_calc$CV=sqrt(hyper_var_all$NME_hypervar_calc$var)/hyper_var_all$NME_hypervar_calc$mean
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="CV")
# 0.1984991
dist_plot_run(as.data.table(hyper_var_all$MML_hypervar_calc),theme_glob,ylab="MML",stat_in="hypervar_logvar",dir=NME_MAV_human_out_dir)
#-0.107
dist_plot_run(as.data.table(hyper_var_all$MML_hypervar_calc),theme_glob,ylab="MML",stat_in="mean",dir=NME_MAV_human_out_dir)
dist_plot_run(as.data.table(hyper_var_all$MML_hypervar_calc),theme_glob,ylab="MML",stat_in="var",dir=NME_MAV_human_out_dir)
hyper_var_all$MML_hypervar_calc$CV=sqrt(hyper_var_all$MML_hypervar_calc$var)/hyper_var_all$MML_hypervar_calc$mean
dist_plot_run(as.data.table(hyper_var_all$MML_hypervar_calc),theme_glob,ylab="MML",stat_in="CV",dir=NME_MAV_human_out_dir)

