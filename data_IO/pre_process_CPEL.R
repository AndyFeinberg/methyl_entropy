rm(list=ls())
source("mainFunctions_sub.R")
setwd('../')
# get all hg19 CpG site ---------------------------------------------------
CpG_hg19=getCpgSitesH19()
saveRDS(CpG_hg19,'../downstream/input/human_analysis/CpG_hg19.rds')
subjects=c("H9","HUES64","skin03","STL001","STL002","STL003",
           "STL011","H1","HuFGM02","112","149","150")
# reading in vcf files ----------------------------------------------------
#Linux for converting vcf to vcf.gz in ../downstream/data/vcfFiles/:
#for fn in *.vcf; do bgzip -c $fn > $fn.gz; done
variant_HetCpG=lapply(subjects,function(x) extractHetCpG('../downstream/data/vcfFiles/',x)) 
names(variant_HetCpG)=subjects
saveRDS(variant_HetCpG,variant_HetCpG_file)

# reading in stat for each sample -------------------------------
#Note here we're using coverage cutoff=5 and boundary check == true
GR=import.subject('../downstream/data/bedGraph_diff/')
saveRDS(GR,GR_file)
GR_allele=import.subject('../downstream/data/bedGraph_allele/',calc='allele')
saveRDS(GR_allele,GR_allele_file)

# reading in analyzed regions for each sample -----------------------------
gff_in=readAllGff('../downstream/data/gff_file/',subjects)
saveRDS(gff_in,gff_in_file)

# Sanity check output regions not in original gff file --------------------
for(subj in subjects){
  cat(paste(subj,':\n'))
  cat(length(subsetByOverlaps(GR_allele[GR_allele$Subject==subj],gff_in[gff_in$Subject==subj],type='equal'))-
        length(GR_allele[GR_allele$Subject==subj]),'\n')
  cat(length(subsetByOverlaps(GR[GR$Subject==subj],gff_in[gff_in$Subject==subj],type='equal'))-length(GR[GR$Subject==subj]),'\n')
}
#Result all 0

# # Extracting genomic features ---------------------------------------------

# genomic_features=getGeneralFeats_CpG("../downstream/input/human_analysis/")
# saveRDS(genomic_features,genomic_features_file)


# creating merged object --------------------------------------------------
GR_merge=GRanges()
for(sp in unique(GR$Sample)){
  cat("Processing",sp,'\n')
  ts=gsub('.*- ','',sp)
  GR_merge=c(GR_merge,stat_merge(GR[GR$Sample==sp],
                                 GR_allele[GR_allele$Sample==sp],
                                 variant_HetCpG[[ts]],CpG_hg19))
}

# Counting splitting events -----------------------------------------------
#Note Jordi splitted regions with large CpG numbers into 2 regions. 
#There's possibility that one region have SNP the other region does not have SNP but we can still separate allele 
#Because allele-separation happened before splitting regions
#The event is rare and only happened in small portion 0.07% of the regions
#They're treated as allele without CG difference in CpG analysis and without binding site difference in motif analysis 
#Per-sample data, highest is 0.56%:   Brain_substantia_nigra_paired - 112
as.data.table(mcols(GR_merge))[,sum(is.na(g1CG))/length(g1CG),by=list(Sample)]
GR_merge$g1CG[is.na(GR_merge$g1CG)]<-GR_merge$g2CG[is.na(GR_merge$g2CG)]<-GR_merge$refCG[is.na(GR_merge$refCG)]<-GR_merge$altCG[is.na(GR_merge$altCG)]<-0
GR_merge=GR_merge[GR_merge$N>=2&!(GR_merge$Sample%in%(c('rep1 - H1','rep2 - H1')))]

# Adding gene information to object ---------------------------------------
GR_merge=add_gene_GR(GR_merge,genomic_features$promoter,'genes_promoter')
GR_merge=add_gene_GR(GR_merge,genomic_features$`gene body`,'genes_body')
GR_merge=add_gene_GR(GR_merge,genomic_features$TSS,'TSS')
GR_merge$Sample[GR_merge$Sample=="merged - H1"] = "ESC - H1"
# loading the varibility file for each sample -----------------------------

agnostic_dir="../downstream/input/human_analysis/NME_expression_var/scRNA/"
GR_merge$hyper_var_fn=NA
GR_merge$tissue[GR_merge$tissue=='Adipose_Tissue_single']='Adipose_single'
GR_merge$hyper_var_fn[GR_merge$Sample %in% c('42_embryonic_stem_cell_single - H9','stem_27_undifferentiated_paired - HUES64',"ESC - H1")]=
  paste0(agnostic_dir,'HESC_1.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Adipose_single']=paste0(agnostic_dir,'AdultAdipose_1.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Bladder_single']=
  paste0(agnostic_dir,'AdultBladder_1.rds',';',agnostic_dir,'AdultBladder_2.rds')
GR_merge$hyper_var_fn[GR_merge$tissue=='Small_Intestine_single']=paste0(agnostic_dir,'AdultIleum_2.rds')
GR_merge$hyper_var_fn[GR_merge$tissue=='Gastric_single']=
  paste0(agnostic_dir,'AdultStomach_2.rds',';',agnostic_dir,'AdultStomach_1.rds',';',
                                  agnostic_dir,'AdultStomach_3.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Left_Ventricle_single']=
  paste0(agnostic_dir,'AdultHeart_1.rds',';',agnostic_dir,'AdultHeart_2.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Lung_single']=
  paste0(agnostic_dir,'AdultLung_1.rds',';',agnostic_dir,'AdultLung_2.rds',';',
                                  agnostic_dir,'AdultLung_3.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Psoas_Muscle_single']=paste0(agnostic_dir,'AdultMuscle_1.rds')
GR_merge$hyper_var_fn[GR_merge$tissue=='Sigmoid_Colon_single']=paste0(agnostic_dir,'AdultColon_1.rds')
GR_merge$hyper_var_fn[GR_merge$tissue=='Spleen_single']=paste0(agnostic_dir,'AdultSpleen_1.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Adrenal_Gland_single']=
  paste0(agnostic_dir,'AdultAdrenalGland_2.rds',';',agnostic_dir,'AdultAdrenalGland_3.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Aorta_single']=paste0(agnostic_dir,'AdultArtery_1.rds')


GR_merge$hyper_var_fn[GR_merge$tissue=='Esophagus_single']=
  paste0(agnostic_dir,'AdultEsophagus_1.rds',';',agnostic_dir,'AdultEsophagus_2.rds')
GR_merge$hyper_var_fn[GR_merge$tissue=='Pancreas_single']=paste0(agnostic_dir,'AdultPancreas_1.rds')

GR_merge$hyper_var_fn[GR_merge$tissue=='Liver_single']=
  paste0(agnostic_dir,'AdultLiver_1.rds',';',agnostic_dir,'AdultLiver_2.rds',';',
                                  agnostic_dir,'AdultLiver_4.rds')

saveRDS(GR_merge,GR_merge_file)

# add CpG information how many CG in each allele etc ----------------------

GR_merge_CpG=GRanges()
for (subj in subjects){GR_merge_CpG=c(GR_merge_CpG,hetCGallele_merged(subj,GR_merge,CpG_hg19,variant_HetCpG,gene_size=500))}
GR_merge_CpG$density=GR_merge_CpG$CG_hg19_extend/GR_merge_CpG$gff_size_extend
GR_merge_CpG$density_diff=(GR_merge_CpG$CG_allele_extend_g1-GR_merge_CpG$CG_allele_extend_g2)/
  GR_merge_CpG$gff_size_extend
saveRDS(GR_merge_CpG,GR_merge_file)

# Processing variant based result -----------------------------------------

variant_HetCpG_meta=fastDoCall('c',lapply(names(variant_HetCpG),variant_meta,variant_in=variant_HetCpG,GR_in=GR_merge))
#Trinucleotide analysis
#variant_HetCpG_meta$mask_tri=unlist(lapply(variant_HetCpG_meta$REF_tri,mask_tri))
saveRDS(variant_HetCpG_meta,variant_HetCpG_meta_file)

# reading in allele-agnostic analysis -------------------------------------
GR_merge=readRDS(GR_merge_file)
in_dir='../downstream/data/allele_agnostic_20kb/'
#all_regions=import.gff3('../downstream/output/human_20kb_allele_agnostic_250bp.gff')
#mcols(all_regions)=mcols(all_regions)[,c("N","score")]
NME_in=GRanges()
MML_in=GRanges()
for(fn in  dir(in_dir)){
  cat('Reading in',fn,'\n')
  stat_in=toupper(sub('.*_','',sub('.bedGraph','',fn)))
  sample_in=sub('_phased.*','',sub('.bedGraph','',fn))
  subject_in=sub('_.*','',sample_in)
  tissue_in=sub(paste0(subject_in,'_'),'',sample_in)
  sample_in=paste0(tissue_in,' - ',subject_in)
  if(sample_in=="ESC_paired - H1"){sample_in="ESC - H1"}
  fn_in=unique(GR_merge[GR_merge$Sample==sample_in]$hyper_var_fn)
if(length(fn_in)==0){fn_in=NA}

  if(stat_in=="NME"){
  #Remove overlapped regions
  NME_in=c(NME_in,read.agnostic(paste0(in_dir,fn),GR_merge[GR_merge$dMML_pval<=pval_cutoff],
                                allele_include = F,sample_in=sample_in,hyper_var_file=fn_in))}
  else if(stat_in=="MML"){
  MML_in=c(MML_in,read.agnostic(paste0(in_dir,fn),GR_merge_in=NULL,allele_include=F,
                                sample_in=sample_in,hyper_var_file=fn_in))}else
   {cat("Error stat_in:", stat_in,'\n')}
}

NME_in$NME=NME_in$score
NME_in=NME_in[NME_in$N>=2]
NME_in=NME_in[!(NME_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
#number of regions check:77457572
saveRDS(NME_in,NME_agnostic_file)
MML_in$MML=MML_in$score
MML_in=MML_in[!(MML_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
MML_in=MML_in[MML_in$N>=2]
#number check: 77460740
saveRDS(MML_in,MML_agnostic_file)

# DNase and control ------------------------------------------------------------
GR_merge=readRDS(GR_merge_file)
in_dir='../downstream/data/agnostic_DNase/'
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
  
  fn_in=unique(GR_merge[GR_merge$Sample==sample_in]$hyper_var_fn)
  if(length(fn_in)==0){fn_in=NA}
  if(stat_in=="NME"&sum(GR_merge$Sample==sample_in)>0){
    NME_in=c(NME_in,read.agnostic(paste0(in_dir,fn),GR_merge[GR_merge$dMML_pval<=pval_cutoff],
                                  allele_include = F,sample_in=sample_in,hyper_var_file=fn_in))}
  else if(stat_in=="MML"&sum(GR_merge$Sample==sample_in)>0){
    MML_in=c(MML_in,read.agnostic(paste0(in_dir,fn),GR_merge_in=NULL,allele_include = F,sample_in=sample_in,hyper_var_file=fn_in))}else
    {cat("Error stat_in:", stat_in,'in sample:',sample_in,'\n')}
}

NME_in$NME=NME_in$score
MML_in$MML=MML_in$score
NME_in=NME_in[!(NME_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
MML_in=MML_in[!(MML_in$Sample%in%c("rep1 - H1","rep2 - H1"))]
NME_in=NME_in[NME_in$N>=2]
MML_in=MML_in[MML_in$N>=2]
#check: NME 57678420
#check: MML 57681712
saveRDS(NME_in,NME_agnostic_DNase_file)
saveRDS(MML_in,MML_agnostic_DNase_file)

# allele-specific regions with agnostic --------------------------------------------------------
in_dir='../downstream/data/allele_specific_region_agnostic/'
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
    NME_in_sp=read.bedGraph.informME(paste0(in_dir,fn))
    NME_in_sp$Sample=sample_in
    NME_in_sp$statistics=stat_in
    NME_in=c(NME_in,NME_in_sp)}
  else if(stat_in=="MML"){
    MML_in_sp=read.bedGraph.informME(paste0(in_dir,fn))
    MML_in_sp$Sample=sample_in
    MML_in_sp$statistics=stat_in
    MML_in=c(MML_in,MML_in_sp)}else
    {cat("Error stat_in:", stat_in,'\n')}
}
rm(NME_in_sp)
rm(MML_in_sp)
#NME in check: 51459969
#MML in check: 51459969
saveRDS(NME_in[NME_in$N>=2],NME_agnostic_ASM_file)
saveRDS(MML_in[MML_in$N>=2],MML_agnostic_ASM_file)

# Allele-agnostic analysis for rest of regions --------------------------------------------------------
in_dir='../downstream/data/compliment_MML_NME_human/'
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
    NME_in_sp=read.bedGraph.informME(paste0(in_dir,fn))
    NME_in_sp$Sample=sample_in
    NME_in_sp$statistics=stat_in
    NME_in=c(NME_in,NME_in_sp)}
  else if(stat_in=="MML"){
    MML_in_sp=read.bedGraph.informME(paste0(in_dir,fn))
    MML_in_sp$Sample=sample_in
    MML_in_sp$statistics=stat_in
    MML_in=c(MML_in,MML_in_sp)}else
    {cat("Error stat_in:", stat_in,'\n')}
}
rm(NME_in_sp)
rm(MML_in_sp)
#NME in check: 202721894
#MML in check: 202721894
saveRDS(NME_in[NME_in$N>=2],NME_agnostic_comp_file)
saveRDS(MML_in[MML_in$N>=2],MML_agnostic_comp_file)


saveRDS(c(readRDS(NME_agnostic_file),
          readRDS(NME_agnostic_DNase_file),
          #readRDS(NME_agnostic_ASM_file),
          readRDS(NME_agnostic_comp_file)
                  ),NME_agnostic_all_file)
saveRDS(c(readRDS(MML_agnostic_file),
          readRDS(MML_agnostic_DNase_file),
          #readRDS(MML_agnostic_ASM_file),
          readRDS(MML_agnostic_comp_file)
),MML_agnostic_all_file)
#Unique analyzed region
unique_gr=unique(granges(MML_all))
CG_hg19=getCpgSitesH19()
CG_hg19_autosome=CG_hg19[seqnames(CG_hg19) %in% c(paste0('chr',c(1:22,'X','Y')))]
length(subsetByOverlaps(CG_hg19_autosome,unique_gr,minoverlap=2))/length(CG_hg19_autosome)#0.7910136
CG_hg19_covered=countOverlaps(CG_hg19_autosome,unique_gr,minoverlap=2)
CG_hg19_covered_tb=table(CG_hg19_covered)
CG_hg19_covered_tb_cov=CG_hg19_covered_tb[names(CG_hg19_covered_tb)!="0"]
CG_hg19_covered_tb_cov[1]/sum(CG_hg19_covered_tb_cov)#0.7925


# reading in mouse MML and NME --------------------------------------------
#Complimentary regions
dir_comp='../downstream/data/compliment_MML_NME_model_mouse/'
MML_in=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*mml"),
                               read.agnostic.mouse,in_dir=dir_comp,mc.cores=20))
MML_in$MML=MML_in$score
MML_in$score=NULL
NME_in=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*nme"),read.agnostic.mouse,in_dir=dir_comp,mc.cores=20))
NME_in$NME=NME_in$score
NME_in$score=NULL
#Analyzed PRC, DNase and control
dir_analyzed='../downstream/data/DNase_control_PRC_MML_NME_model_mouse/'
MML_in_analyzed=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*mml"),
                               read.agnostic.mouse,in_dir=dir_analyzed,mc.cores=20))
MML_in_analyzed$MML=MML_in_analyzed$score
MML_in_analyzed$score=NULL
NME_in_analyzed=fastDoCall('c',mclapply(dir(dir_comp,pattern=".*nme"),read.agnostic.mouse,in_dir=dir_analyzed,mc.cores=20))
NME_in_analyzed$NME=NME_in_analyzed$score
NME_in_analyzed$score=NULL
#Combine two datasets
NME_in=c(NME_in,NME_in_analyzed)
MML_in=c(MML_in,MML_in_analyzed)
#Convert to matrix: note here I didn't filter out the N<=17 since all NME and MML are intersect with UC regions in later analysis, the filtering in only done in UC
NME_in_matrix=agnostic_matrix_conversion(NME_in[NME_in$N>=2])
MML_in_matrix=agnostic_matrix_conversion(MML_in[MML_in$N>=2],'MML')

saveRDS(NME_in_matrix,NME_matrix_file)
saveRDS(MML_in_matrix,MML_matrix_file)
rm(MML_in)
rm(NME_in)
rm(MML_in_matrix)
rm(NME_in_matrix)
gc()
# reading in mouse UC --------------------------------------------
#Complimentary regions
UC_in_dir='../downstream/data/compliment_UC_non_MDS_mouse/'
UC_in=fastDoCall('c',mclapply(dir(UC_in_dir,pattern = '.*uc.bedGraph'),function(x){UC_in=read.agnostic.mouse.uc(paste(UC_in_dir,x,sep=''))
UC_in$UC=UC_in$score
return(UC_in)},mc.cores=20))
UC_in$tissue=sub('-.*','',UC_in$Sample)
UC_in$Sample=sub('.5-.*-E1','.5-E1',UC_in$Sample)
#DNase,control,PRC
UC_in_dir_analyzed='../downstream/data/DNase_control_PRC_non_MDS_mouse/'
UC_in_analyzed=fastDoCall('c',mclapply(dir(UC_in_dir,pattern = paste0(paste(unique(UC_in$tissue),collapse='|'),'.*uc.bedGraph')),
                                       function(x){
  UC_in=read.agnostic.mouse.uc(paste(UC_in_dir_analyzed,x,sep=''))
UC_in$UC=UC_in$score
return(UC_in)},mc.cores=20))
UC_in_analyzed$Sample=sub('.5-.*-E1','.5-E1',UC_in_analyzed$Sample)
UC_in_analyzed=UC_in_analyzed[UC_in_analyzed$Sample %in% unique(UC_in$Sample)]
#Merging data
UC_in=c(UC_in_analyzed[UC_in_analyzed$N>=2&UC_in_analyzed$N<=17],UC_in[UC_in$N>=2&UC_in$N<=17])
#Convert to matrix
UC_in_matrix_ls=mclapply(unique(UC_in$tissue),function(x) agnostic_matrix_conversion(UC_in[UC_in$tissue==x],'UC'),mc.cores=20)
names(UC_in_matrix_ls)=unique(UC_in$tissue)
saveRDS(UC_in_matrix_ls,UC_in_matrix_ls_file)


# UC for mouse MDS comparison ---------------------------------------------
gff_in_compliment=import.gff3(mouse_compliment_gff_file)
gff_in_compliment=paste0(seqnames(gff_in_compliment),':',start(gff_in_compliment),'-',end(gff_in_compliment))
UC_in_MDS_comp=data.table(region=gff_in_compliment)
compliment_MDS_dir='../downstream/data/compliment_UC_MDS_mouse/'
UC_in_MDS_comp_UC=fastDoCall('cbind',
                             mclapply(dir(compliment_MDS_dir,pattern = '.*uc.bedGraph'),function(x){
                               read.agnostic.mouse.uc(paste(compliment_MDS_dir,x,sep=''),matrix=T,fileter_N=2,gff_in=gff_in_compliment)},mc.cores=10))

UC_in_MDS_comp=cbind(UC_in_MDS_comp,UC_in_MDS_comp_UC)
#Filter based on N first to save space
DNase_conrol_MDS_dir='../downstream/data/DNase_control_PRC_MDS_mouse/'
gff_in_DNase=import.gff3(mouse_DNase_control_gff_file)
gff_in_DNase=paste0(seqnames(gff_in_DNase),':',start(gff_in_DNase),'-',end(gff_in_DNase))
UC_in_analyzed_MDS=data.table(region=gff_in_DNase)
UC_in_analyzed_MDS_UC=fastDoCall('cbind',
                                 mclapply(dir(compliment_MDS_dir,pattern = '.*uc.bedGraph'),function(x){
                                   read.agnostic.mouse.uc(paste(DNase_conrol_MDS_dir,x,sep=''),matrix=T,fileter_N=2,gff_in=gff_in_DNase)},mc.cores=10))
UC_in_analyzed_MDS=cbind(UC_in_analyzed_MDS,UC_in_analyzed_MDS_UC)
UC_in_MDS_all=rbind(UC_in_MDS_comp,UC_in_analyzed_MDS)
saveRDS(UC_in_MDS_all,UC_in_MDS_all_file)

# UC for mouse MDS comparison with P0---------------------------------------------
#Note old and new run have different names before and after -vs-
#Fixing this issue
for(fn in dir(compliment_MDS_dir_P0)){
    sample_name=unlist(strsplit(gsub('_uc.bedGraph','',fn),'-vs-'))
    sample_name_rev=paste0(DNase_conrol_MDS_dir,sample_name[2],'-vs-',sample_name[1],'_uc.bedGraph')
    if(file.exists(sample_name_rev)){
         cat('Reversing file name:',fn,' to ',sample_name_rev)
         file.rename(sample_name_rev,paste0(DNase_conrol_MDS_dir,fn))

    }
}
gff_in_compliment=import.gff3(mouse_compliment_gff_file)
gff_in_compliment=paste0(seqnames(gff_in_compliment),':',start(gff_in_compliment),'-',end(gff_in_compliment))
UC_in_MDS_comp_P0=data.table(region=gff_in_compliment)

UC_in_MDS_comp_P0_UC=fastDoCall('cbind',
                             mclapply(dir(compliment_MDS_dir_P0,pattern = '.*uc.bedGraph'),function(x){
                               read.agnostic.mouse.uc(paste(compliment_MDS_dir_P0,x,sep=''),matrix=T,fileter_N=2,gff_in=gff_in_compliment)},mc.cores=20))
UC_in_MDS_comp_P0_UC=cbind(UC_in_MDS_comp_P0,UC_in_MDS_comp_P0_UC)


gff_in_DNase=import.gff3(mouse_DNase_control_gff_file)
gff_in_DNase=paste0(seqnames(gff_in_DNase),':',start(gff_in_DNase),'-',end(gff_in_DNase))
UC_in_analyzed_MDS_P0=data.table(region=gff_in_DNase)
UC_in_analyzed_MDS_P0_UC=fastDoCall('cbind',
                                 mclapply(dir(compliment_MDS_dir_P0,pattern = '.*uc.bedGraph'),function(x){
                                   read.agnostic.mouse.uc(paste(DNase_conrol_MDS_dir,x,sep=''),matrix=T,fileter_N=2,gff_in=gff_in_DNase)},mc.cores=20))
UC_in_analyzed_MDS_P0=cbind(UC_in_analyzed_MDS_P0,UC_in_analyzed_MDS_P0_UC)
UC_in_MDS_all_P0=rbind(UC_in_MDS_comp_P0_UC,UC_in_analyzed_MDS_P0)

saveRDS(UC_in_MDS_all_P0,UC_in_MDS_all_P0_file)
UC_in_all=readRDS(UC_in_MDS_all_file)
UC_in_MDS_all_P0_all=cbind(UC_in_MDS_all_P0,UC_in_all[,-1])
saveRDS(UC_in_MDS_all_P0_all, UC_in_MDS_all_P0_all_file)
# created merged object for all UC, dMML and dNME ----------------------------------------

mml <- readRDS(MML_matrix_file)
mml=convert_GR(mml,direction="matrix")
nme <- readRDS(NME_matrix_file)
nme=convert_GR(nme,direction="matrix")
uc=readRDS(UC_in_matrix_ls_file)
uc=lapply(uc,convert_GR,direction="matrix")
UC_merge=lapply(names(uc),function(x){
  uc_in=uc[[x]]
  mml_in=mml[,grepl(x,colnames(mml))]
  nme_in=nme[,grepl(x,colnames(nme))]
  regions=intersect(intersect(rownames(uc_in),rownames(mml_in)),rownames(nme_in))
  uc_in=uc_in[regions,]
  
  colnames(nme_in)=gsub('.*-','',gsub("-all","",colnames(nme_in)))
  colnames(mml_in)=gsub('.*-','',gsub("-all","",colnames(mml_in)))
  mml_in=mml_in[regions,]
  nme_in=nme_in[regions,]
  time_series=gsub(paste0(x,'-'),'',gsub("-all","",colnames(uc_in)))
  dnme=do.call(cbind,lapply(time_series,function(x){
    return(abs(nme_in[,gsub('-.*','',x)]-nme_in[,gsub('.*-','',x)]))
    
  }))
  dmml=do.call(cbind,lapply(time_series,function(x){
    return(abs(mml_in[,gsub('-.*','',x)]-mml_in[,gsub('.*-','',x)]))
    
  }))
  colnames(uc_in)=paste0("UC-",colnames(uc_in))
  colnames(dmml)=paste0("dMML-",time_series)
  colnames(dnme)=paste0("dNME-",time_series)
  return(cbind(uc_in,dmml,dnme))
})
names(UC_merge)=names(uc)
saveRDS(UC_merge,UC_merge_file)
UC_merge_max_loc=lapply(UC_merge,function(x){
  cat("Percent all data:",sum(rowSums(is.na(x))==0)/nrow(x),'\n')
  x=as.data.frame(x[rowSums(is.na(x))==0,])
  uc_dt=  x[,grepl("UC-",colnames(x))]
  dNME_dt=  x[,grepl("dNME-",colnames(x))]
  dMML_dt=  x[,grepl("dMML-",colnames(x))]
  
  x$dMML_max_pair=apply(dMML_dt,1,max)
  x$dNME_max_pair=apply(dNME_dt,1,max)
  x$UC_max_pair=apply(uc_dt,1,max)
  x$dMML_max_time=gsub('dMML-','',colnames(dMML_dt)[apply(dMML_dt,1,which.max)])
  x$dNME_max_time=gsub('dNME-','',colnames(dNME_dt)[apply(dNME_dt,1,which.max)])
  
  uc_max=apply(uc_dt,1,which.max)
  x$UC_max_time=gsub('UC-','',colnames(uc_dt)[uc_max])
  x$dNME_max_UC_pair=dNME_dt[cbind(seq_along(uc_max), uc_max)]
  #x$UC_max_UC_pair=uc_dt[cbind(seq_along(uc_max), uc_max)]
  x$dMML_max_UC_pair=dMML_dt[cbind(seq_along(uc_max), uc_max)]
  adj_time=paste0(paste0("E",10:15,'.5'),'-',paste0("E",11:16,'.5'))
  uc_max_adj=unlist(apply(x[,(grep(paste0('UC-.*',adj_time,collapse="|",sep=''),colnames(x)))],1,which.max))
  
  x$UC_max_time_adj=gsub('UC-','',colnames(x))[(grepl(paste0('UC-.*',adj_time,collapse="|",sep=''),colnames(x)))][uc_max_adj]
  x$dNME_max_UC_pair_adj=x[,(grepl(paste0('dNME-.*',adj_time,collapse="|",sep=''),colnames(x)))][cbind(seq_along(uc_max_adj), uc_max_adj)]
  x$UC_max_UC_pair_adj=x[,(grepl(paste0('UC-.*',adj_time,collapse="|",sep=''),colnames(x)))][cbind(seq_along(uc_max_adj), uc_max_adj)]
  x$dMML_max_UC_pair_adj=x[,(grepl(paste0('dMML-.*',adj_time,collapse="|",sep=''),colnames(x)))][cbind(seq_along(uc_max_adj), uc_max_adj)]
  return(x)
  
})
names(UC_merge_max_loc)=names(UC_merge)
saveRDS(UC_merge_max_loc,UC_merge_max_loc_file)
#Filtering complete data and regions having all data: should be 4876367
UC_merge=readRDS(UC_merge_file)
UC_merge_max_loc=readRDS(UC_merge_max_loc_file)
UC_merge=UC_filtering(UC_merge)
UC_merge_max_loc=UC_merge=UC_filtering(UC_merge_max_loc)
saveRDS(UC_merge_max_loc,UC_merge_max_loc_file)
saveRDS(UC_merge,UC_merge_file)
#Subset regions by UC>0.1
cluster=readRDS(paste0(dir_cluster_in_01,'uc_0.1_1.rds'))
UC_merge_max_loc_sub=lapply(names(UC_merge_max_loc),function(x) {
  print(x)
  return(UC_merge_max_loc[[x]][names(cluster[[x]]),])
  
})
names(UC_merge_max_loc_sub)=names(UC_merge_max_loc)
saveRDS(UC_merge_max_loc_sub,UC_merge_max_loc_01_file)

#Read in mouse NME and scRNA
NME_in=readRDS(NME_matrix_file)
#From JASON

mcols(NME_in)=mcols(NME_in)[,grepl('limb',colnames(mcols(NME_in)))]

gtf <- fread('../downstream/input/mouse_analysis/grcm38.gtf',data.table = F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
genes <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
genes$gene_name <- gn
NME_in=dist_calc(NME_in,genes)
#Percent gene covered?
length(unique(NME_in[abs(NME_in$dist)<=3000]$gene))/length(genes[seqnames(genes)!="chrM"])#96%
NME_in_dt=convert_GR(NME_in,dir='DT')
NME_in_dt=melt.data.table(NME_in_dt,id.var=c('dist','gene','region'),value.name='NME',variable.name='stage')

NME_in_dt$hyper_var=-100
NME_in_dt$var=-100
NME_in_dt$mean=-100
for(st in unique(NME_in_dt$stage)){
  tt1=proc.time()[[3]]
  if(file.exists(paste0(dir_scRNA_mouse,gsub('E|limb\\.|\\.all','',st),'.rds'))){
    scRNA_in=readRDS(paste0(dir_scRNA_mouse,gsub('E|limb\\.|\\.all','',st),'.rds'))
    scRNA_in=scRNA_in[rownames(scRNA_in)%in% unique(c(NME_in_dt[(stage==st)]$gene)),]
    if(nrow(scRNA_in)>0){
      #Add hypervar to TSS 
      NME_in_dt[(stage==st)]$hyper_var=scRNA_in[NME_in_dt[(stage==st)]$gene,"hypervar_logvar"]
      NME_in_dt[(stage==st)]$var=scRNA_in[NME_in_dt[(stage==st)]$gene,"var"]
      NME_in_dt[(stage==st)]$mean=scRNA_in[NME_in_dt[(stage==st)]$gene,"mean"]
      
    }
  }else{cat("File not exist for ",st,'\n')}
  cat('Finish processing ',sub('E','',st),'in: ',proc.time()[[3]]-tt1,'\n')
  
}
saveRDS(NME_in_dt,NME_mouse_MAV_fn)




