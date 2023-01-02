source('mainFunctions_sub.R')
in_dir='../downstream/input/human_analysis/pseComb/'
NME_in=GRanges()
MML_in=GRanges()
for(fn in  dir(in_dir,pattern="[mn]m[le].bedGraph")){
  cat('Reading in',fn,'\n')
  stat_in=toupper(sub('.*_','',sub('.bedGraph','',fn)))
  sample_in=sub('_phased.*','',sub('.bedGraph','',fn))
  subject_in=sub('_.*','',sample_in)
  tissue_in=sub(paste0(subject_in,'_'),'',sample_in)
  sample_in=paste0(tissue_in,' - ',subject_in)
  
  if(sample_in=="ESC_paired - H1"){sample_in="ESC - H1"}
  
  
  if(stat_in=="NME"){
  NME_in_sp=read.agnostic(paste0(in_dir,fn),GR_merge_in=NULL,
                                allele_include = F,sample_in=sample_in,hyper_var_file=NA)
  NME_in_sp$tissue=tissue_in
  NME_in=c(NME_in,NME_in_sp)
 }
  else if(stat_in=="MML"){
  MML_in_sp=read.agnostic(paste0(in_dir,fn),GR_merge_in=NULL,allele_include = F,sample_in=sample_in,hyper_var_file=NA)
  MML_in_sp$tissue=tissue_in
  MML_in=c(MML_in,MML_in_sp)}else
    {cat("Error stat_in:", stat_in,'\n')}
}
rm(NME_in_sp)
rm(MML_in_sp)
#NME in check: 202721894
#MML in check: 202721894
saveRDS(convert_GR(NME_in[NME_in$N>=2],direction="DT"),"../downstream/output/human_analysis/CPEL_outputs/NME_agnostic_pse.rds")

cor_pse_acutal_plot<-function(statType){
  theme_plot=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=24),
                 axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                 axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                 axis.text.x=element_text(size=16),
                 axis.text.y=element_text(size=16))+theme(legend.position="bottom")
  #Looking at the reverse
  dat_in=readRDS(paste0("../downstream/output/human_analysis/CPEL_outputs/",statType,"_agnostic_pse.rds"))
  dat_all=readRDS(paste0("../downstream/output/human_analysis/CPEL_outputs/",statType,"_agnostic_all.rds"))
  mix_dt=data.table(tissueName=
                  c("endoerm_01_mesoderm_09_paired","endoerm_03_mesoderm_07_paired",
                      "endoerm_05_mesoderm_05_paired", "endoerm_07_mesoderm_03_paired", "endoerm_09_mesoderm_01_paired"),
                    endoermProp=c(0.1,0.3,0.5,0.7,0.9),
                    mesodermProp=c(0.9,0.7,0.5,0.3,0.1)

  )
  dat_in=cbind(dat_in,mix_dt[match(dat_in$tissue,tissueName)])
  dat_all=convert_GR(dat_all[dat_all$Sample %in% c("endoerm_27_paired - HUES64","mesoderm_23_paired - HUES64")],direction="DT")
  dat_in$pureMesoderm=dat_all[Sample=="mesoderm_23_paired - HUES64"][match(dat_in$region,region)]$score
  dat_in$pureEndoerm=dat_all[Sample=="endoerm_27_paired - HUES64"][match(dat_in$region,region)]$score
  dat_in=dat_in[!is.na(pureMesoderm)&!is.na(pureEndoerm)]
  dat_in$predictedEndoerm=(dat_in$score-dat_in$pureMesoderm*dat_in$mesodermProp)/dat_in$endoermProp
  dat_in$predictedMesoderm=(dat_in$score-dat_in$pureEndoerm*dat_in$endoermProp)/dat_in$mesodermProp
  dat_in$predicteMix=dat_in$pureMesoderm*dat_in$mesodermProp+dat_in$pureEndoerm*dat_in$endoermProp

  pdf(paste0("../downstream/output/human_analysis/cell_deconvolution/",statType,"_mesoderm_pred.pdf"),width=5,height=6)
  print(ggplot(dat_in,aes(x=predictedMesoderm,y=pureMesoderm))+geom_bin2d(bins=100)+xlab(paste0("predicted ",statType))+ylab(paste0("actual ",statType))+
      ggtitle(paste0("correlation:",round(cor(dat_in$predictedMesoderm,dat_in$pureMesoderm),digits=3)))+xlim(c(0,1))+scale_fill_gradient(trans="log10")+theme_plot)
  dev.off()
  pdf(paste0("../downstream/output/human_analysis/cell_deconvolution/",statType,"_endoerm_pred.pdf"),width=5,height=6)
  print(ggplot(dat_in,aes(x=predictedEndoerm,y=pureEndoerm))+geom_bin2d(bins=100)+xlab(paste0("predicted ",statType))+ylab(paste0("actual ",statType))+
      ggtitle(paste0("correlation:",round(cor(dat_in$predictedEndoerm,dat_in$pureEndoerm),digits=3)))+xlim(c(0,1))+scale_fill_gradient(trans="log10")+theme_plot)
  dev.off()
  pdf(paste0("../downstream/output/human_analysis/cell_deconvolution/",statType,"_mix_pred.pdf"),width=5,height=6)
  print(ggplot(dat_in,aes(x=predicteMix,y=score))+geom_bin2d(bins=100)+xlab(paste0("predicted ",statType))+ylab(paste0("actual ",statType))+
      ggtitle(paste0("correlation:",round(cor(dat_in$predicteMix,dat_in$score),digits=3)))+xlim(c(0,1))+scale_fill_gradient(trans="log10")+theme_plot)
  dev.off()
  return(dat_in)
}
NME_out=cor_pse_acutal_plot("NME")
MML_out=cor_pse_acutal_plot("MML")
#Id large dNME regions

#Check correlation between psedu NME and MML
statType="NME"
dat_all=readRDS(paste0("../downstream/output/human_analysis/CPEL_outputs/",statType,"_agnostic_all.rds"))
dat_all=convert_GR(dat_all[dat_all$Sample %in% c("endoerm_27_paired - HUES64","mesoderm_23_paired - HUES64")],direction="DT")
dat_all=dcast.data.table(dat_all,region~Sample,value.var="score",fun.aggregate="mean")