# Ken motif processing ----------------------------------------------------
rm(list=ls())
source('mainFunctions_sub.R')
# DNase vs control agnostic --------------------------------------------------------
DNase=readRDS(DNase_mm10_file)
control=readRDS(control_mm10_file)
NME_in=readRDS(NME_matrix_file)
JASPAR_motif=readRDS(JASPAR_motif_mm10_file)
melt_mm10<-function(GR_in,stat_in='NME'){
    GR_in_dt=convert_GR(GR_in,direction='DT')
    GR_in_dt = melt.data.table(GR_in_dt,id.var='region',value.name=stat_in,variable.name='Sample')
    GR_in_gr=convert_GR(GR_in_dt$region,direction='GR')
    mcols(GR_in_gr)=GR_in_dt
    GR_in_gr$region=NULL
    GR_in_gr$Sample=gsub('\\.','_',GR_in_dt$Sample)
    GR_in_gr=GR_in_gr[!is.na(mcols(GR_in_gr)[[stat_in]])]
    return(GR_in_gr)
}
#DNase NME
NME_in_DNase=subsetByOverlaps(NME_in,DNase,type='equal')
NME_in_DNase=melt_mm10(NME_in_DNase)
NME_in_control=subsetByOverlaps(NME_in,control,type='equal')
NME_in_control=melt_mm10(NME_in_control)
split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in_DNase,stat_in='NME',mc.cores=24)
  saveRDS(JASPAR_motif_sp,paste0(Ken_motif_out_dir,'JASPAR_motif_mm10_NME_',i,'_agnostic_DNase.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
#Control NME
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in_control,stat_in='NME',mc.cores=24)
  saveRDS(JASPAR_motif_sp,paste0(Ken_motif_out_dir,'JASPAR_motif_mm10_NME_',i,'_agnostic_conrol.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}


#DNase MML
MML_in=readRDS(MML_agnostic_DNase_file)
MML_in_DNase=subsetByOverlaps(MML_in,DNase,type='equal')
MML_in_control=subsetByOverlaps(MML_in,control,type='equal')

split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=lapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=MML_in_DNase,stat_in='MML')
  saveRDS(JASPAR_motif_sp,paste0(Ken_motif_out_dir,'JASPAR_motif_mm10_MML_',i,'_agnostic_DNase.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
#control MML
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=MML_in_control,stat_in='MML')
  saveRDS(JASPAR_motif_sp,paste0(Ken_motif_out_dir,'JASPAR_motif_mm10_MML_',i,'_agnostic_control.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
