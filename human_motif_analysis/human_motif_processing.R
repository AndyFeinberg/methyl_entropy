# Ken motif processing ----------------------------------------------------
rm(list=ls())
source('mainFunctions_sub.R')
# DNase vs control agnostic --------------------------------------------------------
DNase=readRDS(DNase_hg19_file)
control=readRDS(control_hg19_file)
JASPAR_motif=readRDS(JASPAR_motif_hg19_file)
NME_in=readRDS(NME_agnostic_DNase_file)


#DNase NME
NME_in_DNase=subsetByOverlaps(NME_in,DNase,type='equal')
NME_in_control=subsetByOverlaps(NME_in,control,type='equal')
split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in_DNase,stat_in='NME',mc.cores=24)
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_analysis/Ken_motif/homogeneous/JASPAR_motif_hg19_NME_',i,'_agnostic_DNase.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
#Control NME
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=NME_in_control,stat_in='NME',mc.cores=24)
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_analysis/Ken_motif/homogeneous/JASPAR_motif_hg19_NME_',i,'_agnostic_Control.rds'))
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
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=MML_in_DNase,stat_in='MML',mc.cores=24)
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_analysis/Ken_motif/homogeneous/JASPAR_motif_hg19_MML_',i,'_agnostic_DNase.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
#control MML
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=MML_in_control,stat_in='MML')
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_analysis/Ken_motif/homogeneous/JASPAR_motif_hg19_MML_',i,'_agnostic_Control.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}


# Use allele-specific analysis (use large memory node 1TB) -------------------------------
GR_merge=readRDS(GR_merge_file)
DNase=readRDS(DNase_hg19_file)
control=readRDS(control_hg19_file)
#This is from Ken
JASPAR_motif=readRDS(JASPAR_motif_hg19_file)

#DNase
GR_merge_DNase=subsetByOverlaps(GR_merge,DNase)
GR_merge_control=subsetByOverlaps(GR_merge,control)

#mean NME
split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=GR_merge_DNase,stat_in='NME',mc.cores=12)
  cat('Motif without 49 columns:', which(unlist(lapply(JASPAR_motif_sp,function(x) ncol(mcols(x))!=49))),'\n')
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_analysis/Ken_motif/allelic_motif_hg19/JASPAR_motif_hg19_DNase_allelic_NME_',i,'.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=mclapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=GR_merge_control,stat_in='NME',mc.cores=12)
  cat('Motif without 49 columns:', which(unlist(lapply(JASPAR_motif_sp,function(x) ncol(mcols(x))!=49))),'\n')
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_analysis/Ken_motif/allelic_motif_hg19/JASPAR_motif_hg19_control_allelic_NME_',i,'.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
#mean MML
split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=lapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=GR_merge_DNase,stat_in='MML')
  cat('Motif without 49 columns:', which(unlist(lapply(JASPAR_motif_sp,function(x) ncol(mcols(x))!=49))),'\n')
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_analysis/Ken_motif/allelic_motif_hg19/JASPAR_motif_hg19_DNase_allelic_MML_',i,'.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}
split_data=cut(1:length(JASPAR_motif),breaks=3,label=FALSE)
for(i in 1:3){
  tt1=proc.time()[[3]]
  JASPAR_motif_sp=lapply(JASPAR_motif[split_data==i],NME_dNME_ken,GR_in=GR_merge_control,stat_in='MML')
  cat('Motif without 49 columns:', which(unlist(lapply(JASPAR_motif_sp,function(x) ncol(mcols(x))!=49))),'\n')
  saveRDS(JASPAR_motif_sp,paste0('../downstream/output/human_analysis/Ken_motif/allelic_motif_hg19/JASPAR_motif_hg19_control_allelic_MML_',i,'.rds'))
  JASPAR_motif_sp=NA
  proc.time()[[3]]-tt1
}