rm(list=ls())
source('mainFunctions_sub.R')
meta_data_sample=as.data.table(read_xlsx('../downstream/input/mouse_analysis/Ecker_TableS1.xlsx',skip=2))
meta_data_sample=meta_data_sample[,list(tissue,stage)]
meta_data_sample$tissue=tolower(meta_data_sample$tissue)
meta_data_sample[tissue=='craniofacial']$tissue='EFP'
meta_data_sample$stage=gsub('E','day',meta_data_sample$stage)
meta_data_sample$stage=gsub('\\.5','_5',meta_data_sample$stage)
meta_data_sample$stage=gsub('P0','day0',meta_data_sample$stage)
tissue_to_analyze=c('forebrain','heart','hindbrain','EFP','limb','liver','midbrain')
meta_data_sample=meta_data_sample[tissue%in%tissue_to_analyze& grepl('day',stage)]
CPEL_sub_dir='../downstream/input/mouse_analysis/CPEL_submission_script/'
#writing submission for model estimation
#submission command:./sbatching_mouse_all_st
model_est =file(paste0(CPEL_sub_dir,'sbatching_mouse_all_st'))

writeLines(c('#!/usr/bin/env bash',
             paste0('sbatch cpelasm_allele_agnostic.slrm mm10 ',meta_data_sample$tissue,'_',meta_data_sample$stage,'_all false')), model_est)
close(model_est)

#writing submission for tissue-specific UC estimation
#submission command sbatch cpelasm_allele_agnostic_uc_array_tissue.slrm UC_submission_all_non_MDS_arg: note to change #SBATCH --array=1-906 number to number of rows in arg file 
meta_data_sample$sample_n=paste0('mm10_',meta_data_sample$tissue,'_',meta_data_sample$stage,'_all')
tissue_submission_arg=c()
for(ts in unique(meta_data_sample$tissue)){
  tissue_submission_arg=c(tissue_submission_arg,
                       apply(combn(meta_data_sample[tissue==ts]$sample_n,2),2,function(x) paste0(x[1],' ',x[2])))
 
    
    
}
UC_non_MDS =file(paste0(CPEL_sub_dir,'UC_submission_all_non_MDS_arg'))

writeLines(sort(tissue_submission_arg), UC_non_MDS)
close(UC_non_MDS)

#writing submission for all MDS estimation: should be 1275 for all

#submission command sbatch cpelasm_allele_agnostic_uc_array.slrm UC_submission_all_MDS_arg_compliment: note to change #SBATCH --array=1-906 number to number of rows in arg file 
#use linux to exclude the ones already in non_MDS ones: comm -23 UC_submission_all_MDS_arg UC_submission_all_non_MDS_arg >UC_submission_all_MDS_arg_non_tissue
UC_MDS =file(paste0(CPEL_sub_dir,'UC_submission_all_MDS_arg'))
writeLines(sort(apply(combn(meta_data_sample$sample_n,2),2,function(x) paste0(x[1],' ',x[2]))), UC_MDS)
close(UC_MDS)
#Day0 only: 34
#grep day0 UC_submission_all_non_MDS_arg >UC_submission_all_non_MDS_arg_day0
#sbatch cpelasm_allele_agnostic_uc_array.slrm UC_submission_all_non_MDS_arg_day0


