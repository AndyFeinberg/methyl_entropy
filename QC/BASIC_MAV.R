 source('mainFunctions_sub.R')
library(BASiCS)
#https://www.bioconductor.org/packages/release/bioc/vignettes/BASiCS/inst/doc/BASiCS.html
#Check for install: https://github.com/stan-dev/rstan/issues/892
Data <- makeExampleBASiCS_Data()
Chain <- BASiCS_MCMC(
  Data = Data,
  N = 100, Thin = 20, Burn = 40,
  PrintProgress = FALSE, Regression = TRUE,WithSpikes = FALSE
)
#WE USE A SMALL NUMBER OF ITERATIONS FOR ILLUSTRATION PURPOSES ONLY. LARGER NUMBER OF ITERATIONS ARE USUALLY REQUIRED TO ACHIEVE CONVERGENCE. OUR RECOMMENDED SETTING IS N=20000, Thin=20 and Burn=10000
HCL_dir='../downstream/data/HCL_scRNA/count/'
fn=dir(HCL_dir)[1]
tt1=proc.time()[[3]]
expr=readRDS(paste0(HCL_dir,fn))
sce <- SingleCellExperiment(list(counts=expr),
    colData=DataFrame(label=colnames(expr)),
    rowData=DataFrame(length=rownames(expr)),
    metadata=list(study=fn)
)
Chain <- BASiCS_MCMC(
  Data = sce,
  N = 1000, Thin = 10, Burn = 50,
  PrintProgress = TRUE, Regression =TRUE,WithSpikes = TRUE
)
chain_summary=Summary(Chain)
res=data.table(
  mean=displaySummaryBASiCS(chain_summary, Param = "mu")[,"median"],
  var=displaySummaryBASiCS(chain_summary, Param = "delta")[,"median"],
  hypervar_logvar=displaySummaryBASiCS(chain_summary, Param = "epsilon")[,"median"],
  gene=names(displaySummaryBASiCS(chain_summary, Param = "epsilon")[,"median"])
)
#6 hours per sample
saveRDS(res,"../downstream/output/human_analysis/QC/BASIC_MAV_N20000_tin20_burn10000.rds")
#Compare data 
res_low_N=readRDS("../downstream/output/human_analysis/QC/BASIC_MAV_N100_tin20_burn40.rds")
res_high_N=readRDS("../downstream/output/human_analysis/QC/BASIC_MAV_N20000_tin20_burn10000.rds")
res_high_N=res_high_N[match(res_low_N$gene,gene)]
res_MAV=readRDS("../downstream/input/human_analysis/NME_expression_var/scRNA/AdultAdipose_1.rds")
cor.test(res_low_N$hypervar_logvar,res_MAV[res_low_N$gene,"hypervar_logvar"])#0.4034069
cor.test(res_high_N$hypervar_logvar,res_MAV[res_high_N$gene,"hypervar_logvar"])#0.5900582
cor.test(res_low_N$hypervar_logvar,res_high_N$hypervar_logvar)
pdf("../downstream/output/human_analysis/QC/BASIC_MAV_N.pdf")
hyperVarDt=data.table(lowN=res_low_N$hypervar_logvar,highN=res_high_N[match(res_low_N$gene,gene)]$hypervar_logvar,MAV=res_MAV[res_high_N$gene,"hypervar_logvar"])
print(ggplot(hyperVarDt,aes(x=highN,y=lowN))+ggtitle("high N vs low N")+geom_smooth()+geom_point(alpha=0.1))
print(ggplot(hyperVarDt,aes(x=MAV,y=lowN))+ggtitle("MAV vs low N")+geom_smooth()+geom_point(alpha=0.1))
print(ggplot(hyperVarDt,aes(x=MAV,y=highN))+ggtitle("MAV vs high N")+geom_smooth()+geom_point(alpha=0.1))
dev.off()



pdf("../downstream/output/human_analysis/QC/BASIC_MAV.pdf")
res=res[!is.na(hypervar_logvar)]
print(ggplot(res,aes(x=log2(mean),y=log2(var)))+ggtitle("var")+geom_smooth()+geom_point(alpha=0.1))
print(ggplot(res,aes(x=log2(mean),y=hypervar_logvar))+ggtitle("hypervar")+geom_smooth()+geom_point(alpha=0.1))
dev.off()
cat("Run finish in ",proc.time()[[3]]-tt1,"\n")
#5 hours per sample?

 
 if(use_spikeIn){
    Tech=grepl("^ERCC",rownames(expr),ignore.case=T)
    print(which(Tech))
    expr_gene=expr[!Tech,]
    spikeIn=expr[Tech,]
    sce <- SingleCellExperiment(list(counts=expr_gene),
        colData=DataFrame(label=colnames(expr_gene)),
        rowData=DataFrame(length=rownames(expr_gene)),
        metadata=list(study=fn)
    )
    altExp(sce,"ERCC")= SingleCellExperiment(list(counts=spikeIn+0.0001),
        colData=DataFrame(label=colnames(spikeIn)),
        rowData=DataFrame(length=rownames(spikeIn)),
        metadata=list(study=fn)
    )
    rowData(altExp(sce))$concentration=2
    Chain <- BASiCS_MCMC(
      Data = sce,
      N = N, Thin = 20, Burn = Burn,
      PrintProgress = TRUE, Regression =TRUE,WithSpikes = TRUE
    )
  }else{


    
saveRDS(res,"../downstream/output/human_analysis/QC/BASIC_MAV_N10_tin20_burn40_spike.rds")
saveRDS(res,"../downstream/output/human_analysis/QC/BASIC_MAV_N10_tin20_burn40_noSpike.rds")
res_MAV=readRDS("../downstream/input/human_analysis/NME_expression_var/scRNA/HESC_1.rds")
res_spike=readRDS("../downstream/output/human_analysis/QC/BASIC_MAV_N10_tin20_burn40_spike.rds")
res_noSpike_large=readRDS("../downstream/input/human_analysis/HCL_scRNA_BASiC/HESC_1.rds")
genes=intersect(rownames(res_MAV),rownames(res))
cor.test(res_MAV[genes,"hypervar_logvar"],res[genes,"hypervar_logvar"])#0.308
genes=intersect(rownames(res_noSpike_large),rownames(res))
cor.test(res_noSpike_large[genes,"hypervar_logvar"],res[genes,"hypervar_logvar"])#0.2814
genes=intersect(rownames(res_noSpike_large),rownames(res_spike))
cor.test(res_noSpike_large[genes,"hypervar_logvar"],res_spike[genes,"hypervar_logvar"])#0.507
genes=intersect(rownames(res_noSpike_large),rownames(res_MAV))
cor.test(res_noSpike_large[genes,"hypervar_logvar"],res_MAV[genes,"hypervar_logvar"])#0.507


#with replicates: adipose tissue, 3
cor.test(MAV[genes,"hypervar_logvar"],res[,"hypervar_logvar"])#0.5375332 
#without replicates: adipose tissue, 3
cor.test(MAV[genes,"hypervar_logvar"],res[,"hypervar_logvar"])#0.5536888 

#This is without combine replicates
hyper_var_dir="../downstream/input/human_analysis/HCL_scRNA_BASiC/"
scRNA_dir="../downstream/input/human_analysis/NME_expression_var/scRNA/"
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
  hyper_var_file=gsub(scRNA_dir,hyper_var_dir,hyper_var_file)
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
NME_hypervar_calc_no_NA=NME_hypervar_calc[!is.na(NME_hypervar_calc$hypervar_logvar)]
NME_hypervar_calc_no_NA_dt=convert_GR(NME_hypervar_calc_no_NA,"DT")
NME_hypervar_calc_no_NA_dt_tss=NME_hypervar_calc_no_NA_dt[abs(dist)<=250]
#NME: 53476752 check
#MML: check: 53478342
saveRDS(list(NME_hypervar_calc=NME_hypervar_calc,
             MML_hypervar_calc=MML_hypervar_calc),
             paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV_BASiC.rds'))
hyper_var_all=readRDS( paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV_BASiC.rds'))
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="hypervar_logvar",dir=NME_MAV_human_out_dir)
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="mean",dir=NME_MAV_human_out_dir)
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="var",dir=NME_MAV_human_out_dir)
hyper_var_all$NME_hypervar_calc$CV=sqrt(hyper_var_all$NME_hypervar_calc$var)/hyper_var_all$NME_hypervar_calc$mean
dist_plot_run(as.data.table(hyper_var_all$NME_hypervar_calc),theme_glob,ylab="NME",stat_in="CV")

hyper_var_all_old=readRDS( paste0(NME_MAV_human_out_dir,'allele_agnostic_var_homogeneous2_MAV.rds'))
hyper_var_all_old_NME=hyper_var_all_old$NME_hypervar_calc
hyper_var_all_old_NME=data.table(Sample=hyper_var_all_old_NME$Sample,gene=hyper_var_all_old_NME$gene,hypervar_logvar=hyper_var_all_old_NME$hypervar_logvar)
hyper_var_all_old_NME= unique(hyper_var_all_old_NME[!is.na(hypervar_logvar)])
hyper_var_all_old_NME$gene_Sample=paste0(hyper_var_all_old_NME$gene,hyper_var_all_old_NME$Sample)

hyper_var_all_NME=hyper_var_all$NME_hypervar_calc
hyper_var_all_NME=data.table(Sample=hyper_var_all_NME$Sample,gene=hyper_var_all_NME$gene,hypervar_logvar=hyper_var_all_NME$hypervar_logvar)
hyper_var_all_NME= unique(hyper_var_all_NME[!is.na(hypervar_logvar)])
hyper_var_all_NME$gene_Sample=paste0(hyper_var_all_NME$gene,hyper_var_all_NME$Sample)
hyper_var_all_NME=hyper_var_all_NME[match(hyper_var_all_old_NME$gene_Sample,gene_Sample)]
cor.test(hyper_var_all_NME$hypervar_logvar,hyper_var_all_old_NME$hypervar_logvar)#0.594