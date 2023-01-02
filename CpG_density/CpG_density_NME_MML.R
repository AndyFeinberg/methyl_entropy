rm(list=ls())
setwd('../')
source("mainFunctions_sub.R")
#Define ggplot theme

# Preprocess SNP files to get unique SNP and trinucleotide -----------------------------------------------------
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
variant_HetCpG_meta_dt=convert_GR(variant_HetCpG_meta,direction='DT')
variant_HetCpG_meta_dt$SNP=apply(variant_HetCpG_meta_dt[,list(REF,ALT)],1,function(x) paste(x,collapse = '->'))
variant_HetCpG_meta_dt$tri_SNP=paste0(variant_HetCpG_meta_dt$REF_tri,'->',variant_HetCpG_meta_dt$ALT_tri)
saveRDS(variant_HetCpG_meta_dt,variant_HetCpG_meta_dt_file)
single_SNP_unique=unique_mutation(unique(variant_HetCpG_meta_dt$SNP))
#manually order unique mutation
single_SNP_unique[single_SNP_unique=="T->C"]="C->T"
single_SNP_unique[single_SNP_unique=="T->G"]="G->T"
single_SNP_unique[single_SNP_unique=="A->G"]="G->A"
single_SNP_unique[single_SNP_unique=="G->C"]="C->G"
single_SNP_unique[single_SNP_unique=="A->T"]="T->A"
single_SNP_unique[single_SNP_unique=="A->C"]="C->A"
mutation_tri_unique=unique_mutation(unique(variant_HetCpG_meta_dt$tri_SNP))
# #correct order based on single SNP mutation
for(i in 1:length(mutation_tri_unique)){
  tri=mutation_tri_unique[i]

  if(!tri_to_SNP(tri)%in%single_SNP_unique){
    tri_1=gsub('->.*','',tri)
    tri_2=gsub('.*->','',tri)
    mutation_tri_unique[i]=paste0(tri_2,'->',tri_1)


  }

}
#Make a hashtable-like to easy access data
mutation_tri_unique_dt=data.table(raw=names(mutation_tri_unique),unique_fw=mutation_tri_unique)
#Merge reverse complement
single_SNP_unique[single_SNP_unique=="G->A"]="C->T"
single_SNP_unique[single_SNP_unique=="G->T"]="C->A"
variant_HetCpG_meta_dt$SNP=single_SNP_unique[variant_HetCpG_meta_dt$SNP]
mutation_tri_unique_dt$rev_comp=unlist(lapply(mutation_tri_unique_dt$unique_fw,reverse_comp_SNP))
mutation_tri_unique_dt$rev_comp_inv=paste0(gsub('.*->','',mutation_tri_unique_dt$rev_comp),'->',
                                           gsub('->.*','',mutation_tri_unique_dt$rev_comp))
mutation_tri_unique_dt$rev_comp_unique="NA"
for(i in 1:nrow(mutation_tri_unique_dt)){
  #Only look for trinucleotide in the single SNP catogry
  tri_in=mutation_tri_unique_dt$unique_fw[i]
  if(tri_to_SNP(tri_in)%in%single_SNP_unique){
    mutation_tri_unique_dt[(rev_comp==tri_in|unique_fw==tri_in|rev_comp_inv==tri_in)&rev_comp_unique=="NA"]$rev_comp_unique=tri_in
    
    
  }
  
  
}
#Should be 52 unique ones
mutation_tri_unique=mutation_tri_unique_dt$rev_comp_unique
names(mutation_tri_unique)=mutation_tri_unique_dt$raw
#Check how many of them have CG in right, should be 0 since the trinucleotide is ordered
gainCG_idx=which(grepl("CG",sub('.*->','',mutation_tri_unique))&!grepl("CG",sub('->.*','',mutation_tri_unique)))#Should only be 0

# calculate the relative NME using the order of SNP -----------------------

variant_HetCpG_meta_dt$tri_SNP_unique=mutation_tri_unique[variant_HetCpG_meta_dt$tri_SNP]
variant_HetCpG_meta_dt$dNME_relative=as.numeric(NA)
#Use minus sign to get NME right-NME left, .e.g (NME in GTG- NME in GCG) for (GCG->GTG)
variant_HetCpG_meta_dt$dNME_relative=-variant_HetCpG_meta_dt[,list(dNME_relative=dNME_relative_calc(genome1_tri,genome2_tri,NME1,NME2,tri_SNP,tri_SNP_unique,SNP)), 
                                                            by = seq_len(nrow(variant_HetCpG_meta_dt))]$dNME_relative
  

#Convert everything to gain CG, id SNPs with CG changes
variant_HetCpG_meta_dt$CpG_change='Lose CpG'
variant_HetCpG_meta_dt[((grepl('CG',REF_tri )) & (grepl('CG',ALT_tri)))|(!grepl('CG',REF_tri )) & (!grepl('CG',ALT_tri))]$CpG_change='No CpG change'
saveRDS(variant_HetCpG_meta_dt,variant_HetCpG_meta_dt_uq_file)
variant_HetCpG_meta_dt=readRDS(variant_HetCpG_meta_dt_uq_file)
variant_HetCpG_meta_dt$CpG_change=gsub('CG','CpG',variant_HetCpG_meta_dt$CpG_change)
#Generating Figure 4B and calculating OR for dNME
#Param initialization & color theme
SNP_all=list()
SNP_het=list()
SNP_box=list()
sig_v=0.75
sig_h_pos=-0.05
sig_h_neg=1.2
color_theme=c(rainbow(length(unique(variant_HetCpG_meta_dt$SNP))))
variant_SNP_tri=data.table()
variant_SNP_tri_out=list()
names(color_theme)=unique(variant_HetCpG_meta_dt$SNP)
theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=12),
                                 axis.title.x=element_text(hjust=0.5,size=9,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=9,face="bold"),
                                 axis.text.x=element_text(size=7),
                                 axis.text.y=element_text(size=7))
# text=element_text(family="Space Mono"))
for (sn in unique(variant_HetCpG_meta_dt$SNP)){
  #OR calculation
  for(tri in unique(variant_HetCpG_meta_dt[SNP==sn]$tri_SNP_unique)){
    variant_SNP_tri_OR=OR_calc(variant_HetCpG_meta_dt[SNP==sn &dNME_pval<=pval_cutoff],tri,"tri_SNP_unique",pval_cutoff)
    variant_SNP_tri_OR$CpG_change=unique(unique(variant_HetCpG_meta_dt[SNP==sn & tri_SNP_unique==tri]$CpG_change))
    
      variant_SNP_tri_OR$CpG_change=unique(unique(variant_HetCpG_meta_dt[SNP==sn & tri_SNP_unique==tri]$CpG_change))
      
      
   
    variant_SNP_tri=rbind(variant_SNP_tri,variant_SNP_tri_OR)
  }
  #get HetCpG
  variant_SNP_tri$CpG_change=factor(variant_SNP_tri$CpG_change,levels=c('Lose CpG','No CpG change'))
  variant_SNP_tri=variant_SNP_tri[order(OR,decreasing=F)]
 
  variant_SNP_tri$SNP=gsub('->','\u2794',variant_SNP_tri$SNP)
  variant_SNP_tri$SNP=factor(variant_SNP_tri$SNP,levels = variant_SNP_tri$SNP)
  
  variant_SNP_tri$FDR=p.adjust(variant_SNP_tri$pvalue,method='BH')
  variant_SNP_tri$significant=add.significance.stars(variant_SNP_tri$FDR, cutoffs = c(0.05, 0.01, 0.001))
  #Plotting
   SNP_het[[sn]]=ggplot(variant_SNP_tri,aes(x=SNP,y=log(OR),fill=CpG_change))+geom_bar(stat="identity")+ylab('')+xlab("")+
    geom_errorbar(aes(ymin=log(lowerCI), ymax=log(upperCI)), width=.4,position=position_dodge(.9),size=0.25)+ggtitle(gsub('->',' \u2794 ',sn))+#ylim(c(0,max(variant_SNP_tri$upperCI)*1.5))+
    theme_glob+theme(legend.position = "bottom",legend.title = element_blank())+
     ylim(c(-1.5,1.5))+
     scale_fill_manual(values=c("No CpG change"="grey","Lose CpG"="light blue"))+
     #geom_text(data=variant_SNP_tri[OR>1],aes(label=significant,y=log(upperCI)*1),vjust =sig_v,hjust=sig_h_pos)+
     #geom_text(data=variant_SNP_tri[OR<1],aes(label=significant,y=log(lowerCI)*1),vjust =sig_v,hjust=sig_h_neg)+
     coord_flip()
   

   variant_SNP_tri_out[[sn]]=variant_SNP_tri
   variant_SNP_tri=data.table()

}
saveRDS(SNP_het,'../downstream/output/human_analysis/CpG_density/SNP_het_CpG.rds')
saveRDS(variant_SNP_tri_out,'../downstream/output/human_analysis/CpG_density/variant_SNP_tri_out_CpG.rds')
#Getting png files with mono-spaced font in windows setting, use the variant_SNP_tri_out
variant_SNP_tri_out=readRDS('../downstream/output/human_analysis/CpG_density/variant_SNP_tri_out_CpG.rds')
library(extrafont)
loadfonts()
SNP_het=readRDS('../downstream/output/human_analysis/CpG_density/SNP_het_CpG.rds')
png('../downstream/output/human_analysis/CpG_density/variant_OR_tri3_two_cat_greater_CG_bg_rev_hg19.png',
    width=7,height=7,units='in',res=1080, family = 'Consolas')
#SNP_het=SNP_het[c("C>G", names(SNP_het)[names(SNP_het)!="C>G"])]
ggarrange(plotlist=lapply(SNP_het,function(x) x+ ylab("log(Odds Ratio)")+theme( axis.title.x=element_text(hjust=0.5,size=16,face="bold"))), 
          nrow=2,ncol=2,common.legend = T,legend="top")
dev.off()

#Calculate OR for SNPs gaining CG, numbers in text
OR_calc(variant_HetCpG_meta_dt[dNME_pval<=pval_cutoff],'Lose CpG',"CpG_change")

# Density analysis using allele-specific way using regions------------------------------
GR_merge=readRDS(GR_merge_file)

GR_merge$CpGdiff=GR_merge$g1CG-GR_merge$g2CG
#CpG density vs dNME,here uses an extended density 
GR_merge$dNME_relative=GR_merge$NME1-GR_merge$NME2
GR_merge$dMML_relative=GR_merge$MML1-GR_merge$MML2
GR_merge_dt=convert_GR(GR_merge,'DT')
#density difference difference
GR_merge_dt$density_diff=GR_merge_dt[,(CG_allele_extend_g1-CG_allele_extend_g2)/CGcont_exp ]
#Correlation:-0.390
cor.test(GR_merge_dt[dNME_pval<=pval_cutoff&GR_merge_dt$CpGdiff!=0]$dNME_relative, 
         GR_merge_dt[dNME_pval<=pval_cutoff&dNME_pval<=pval_cutoff&GR_merge_dt$CpGdiff!=0]$density_diff)

#Figure C for density
#compare allele with more CpG vs allele with fewer CpG
GR_merge_dt_CG_diff=GR_merge_dt[dNME_pval<=pval_cutoff]
GR_merge_dt_CG_diff[,`more CG`:=c(NME1,NME2)[which.max(c(g1CG,g2CG))],by=seq_len(nrow(GR_merge_dt_CG_diff))][,`fewer CG`:=c(NME1,NME2)[which.min(c(g1CG,g2CG))],by=seq_len(nrow(GR_merge_dt_CG_diff))]
GR_merge_dt_CG_diff_mt=melt.data.table(GR_merge_dt_CG_diff[g1CG!=g2CG,list(`more CG`,`fewer CG`)],id.vars=NULL,variable.name = "region_type",value.name="NME")
all_NME=c(GR_merge_dt_CG_diff$NME1,GR_merge_dt_CG_diff$NME2)
GR_merge_dt_CG_diff_mt=rbind(GR_merge_dt_CG_diff_mt,data.table(region_type="all",NME=all_NME))

pdf(paste0(figure_path,'CpG_number_NME_hg19_violin.pdf'),width=4,height=2)
ggplot(GR_merge_dt_CG_diff_mt,aes(x=region_type,y=NME,fill=region_type))+geom_violin()+xlab("allele type")+theme_glob+theme(legend.position = "none")
dev.off()

CG_difference=ggplot(GR_merge_dt_CG_diff[g1CG!=g2CG],aes(x=`more CG`,y=`fewer CG`))+geom_hex()+theme_glob+theme(legend.position = "none")+ggtitle("Allele has CG difference")
noCG_difference=ggplot(GR_merge_dt_CG_diff[g1CG==g2CG],aes(x=NME1,y=NME2))+geom_hex()+theme_glob+theme(legend.position = "none")+ggtitle("Allele has no CG difference")
pdf(paste0(figure_path,'CpG_number_NME_hg19_scatter.pdf'),width=4.75,height=10)
ggarrange(CG_difference,noCG_difference,ncol=1)
dev.off()

t.test(x=GR_merge_dt_CG_diff[g1CG!=g2CG]$`more CG`,y=GR_merge_dt_CG_diff[g1CG!=g2CG]$`fewer CG`,paired=TRUE,alternative="less")
#Categorizing regions
GR_merge_dt$CpG_stat="No difference"
GR_merge_dt[CpGdiff!=0]$CpG_stat="With CpG difference"
GR_merge_dt$CpG_stat=factor(GR_merge_dt$CpG_stat,levels = c("With CpG difference","No difference"))
#Always using allele with more CG minus alleles with less CG
GR_merge_dt$dNME_relative_more_less=GR_merge_dt$dNME_relative
GR_merge_dt[GR_merge_dt$CpGdiff!=0]$dNME_relative_more_less=GR_merge_dt[GR_merge_dt$CpGdiff!=0]$dNME_relative*sign(GR_merge_dt[GR_merge_dt$CpGdiff!=0]$CpGdiff)
#Test for if "With CpG difference" is significantly smaller than "No difference"
t.test(GR_merge_dt[CpGdiff!=0&dNME_pval<=pval_cutoff]$dNME_relative_more_less,alternative="less")
#Figure C
pdf(paste0(figure_path,'CpG_number_NME_hg19.pdf',width=7,height=7))
ggplot(GR_merge_dt[dNME_pval<=pval_cutoff],aes(y=dNME_relative_more_less,x=CpG_stat,fill=CpG_stat))+
  geom_violin()+xlab("")+
  theme_glob+ylab('relative dNME')+theme(legend.position = "none")
dev.off()
#Line plot for dNME vs density
GR_merge_dt_sig_density_diff=GR_merge_dt[dNME_pval<=pval_cutoff&density_diff!=0]
GR_merge_dt_sig_density_diff$density_difference_quantile=ecdf(GR_merge_dt_sig_density_diff$density_diff)(GR_merge_dt_sig_density_diff$density_diff)
pdf(paste0(figure_path,'CpG_density_dNME_ratio_hg19.pdf'),width=7,height=7)
ggplot(GR_merge_dt_sig_density_diff,aes(x=density_difference_quantile,y=dNME_relative))+geom_smooth(fill='light blue')+
  xlab("CpG density ratio quantile")+ylab("relative dNME")+
  theme_glob+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
#Line plot for dMML
GR_merge_dt_sig_density_diff=GR_merge_dt[dMML_pval<=pval_cutoff&density_diff!=0]
GR_merge_dt_sig_density_diff$density_difference_quantile=ecdf(GR_merge_dt_sig_density_diff$density_diff)(GR_merge_dt_sig_density_diff$density_diff)
pdf('../downstream/output/graphs_tables/CpG_density_dMML_ratio_hg19.pdf',width=7,height=7)
ggplot(GR_merge_dt_sig_density_diff,aes(x=density_difference_quantile,y=dMML_relative))+geom_smooth(fill='light blue')+
  xlab("CpG density ratio quantile")+ylab("relative dMML")+
  theme_glob+
  theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#allele-agnostic density ---------------------------------------
CG_exp_agnostic_hg19_file='../downstream/output/human_analysis/CpG_density/analyzed_region_CG_hg19.rds'
#NME
NME_in=readRDS(NME_agnostic_file)
analyzed_region=unique(granges(NME_in))
gr_seq=getSeq(Hsapiens,analyzed_region,as.character=T)
analyzed_region$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
saveRDS(analyzed_region,CG_exp_agnostic_hg19_file)

analyzed_region=readRDS(CG_exp_agnostic_hg19_file)
NME_in=readRDS(NME_agnostic_file)
NME_in_olap=findOverlaps(NME_in,analyzed_region,type='equal')
NME_in$CGcont_exp[queryHits(NME_in_olap)]=analyzed_region$CGcont_exp[subjectHits(NME_in_olap)]
CpG_hg19=readRDS(CpG_hg19_file)
NME_in$CG_hg19=countOverlaps(NME_in,CpG_hg19)
NME_in$density=NME_in$CG_hg19/NME_in$CGcont_exp
cor.test(NME_in$density,NME_in$NME,method='pearson')#-0.21
#Making boxplot of this with different interval
NME_in$density_quant=findInterval(NME_in$density,seq(0,1,0.1))
quant_conv=c(paste0(seq(0,0.9,0.1),'-',seq(0.1,1,0.1)),'>1')
NME_in$density_quant=factor(quant_conv[NME_in$density_quant],levels=quant_conv)
#PLotting boxplot
pdf(paste0(figure_path,'CpG_density_NME_boxplot_CG_exp.pdf'),width=3.5,height=3.5)
ggplot(as.data.frame(mcols(NME_in)),aes(x=density_quant, y=NME))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()   

#Feature enrichment of low NME regions
genomic_features=readRDS(genomic_features_file)
#Figure S5
olap_islands=findOverlaps(NME_in,genomic_features$`CpG island`)
olap_shores=findOverlaps(NME_in,genomic_features$`CpG shore`)
olap_shelf=findOverlaps(NME_in,genomic_features$`CpG shelf`)
olap_open_sea=findOverlaps(NME_in,genomic_features$`CpG open sea`)


CpG_density_NME=rbind(data.table(NME=NME_in$NME[queryHits(olap_islands)],feature='islands'),
                      data.table(NME=NME_in$NME[queryHits(olap_shores)],feature='shores'),
                      data.table(NME=NME_in$NME[queryHits(olap_shelf)],feature='shelf'),
                      data.table(NME=NME_in$NME[queryHits(olap_open_sea)],feature='open sea'))

pdf(paste0(figure_path,'CpG_density_NME_features.pdf'),width=3.5,height=3.5)
ggplot(as.data.frame(CpG_density_NME),aes(x=feature, y=NME))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("NME")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off() 
#MML
MML_in=readRDS(MML_agnostic_file)
CpG_hg19=readRDS(CpG_hg19_file)
MML_in$CG_hg19=countOverlaps(MML_in,CpG_hg19)
analyzed_region=readRDS(CG_exp_agnostic_hg19_file)
MML_in_olap=findOverlaps(MML_in,analyzed_region,type='equal')
MML_in$CGcont_exp[queryHits(MML_in_olap)]=analyzed_region$CGcont_exp[subjectHits(MML_in_olap)]
MML_in$density=MML_in$CG_hg19/MML_in$CGcont_exp
cor.test(MML_in$density,MML_in$MML,method='pearson')#-0.346
#Making boxplot of this with different interval
MML_in$density_quant=findInterval(MML_in$density,seq(0,1,0.1))
quant_conv=c(paste0(seq(0,0.9,0.1),'-',seq(0.1,1,0.1)),'>1')
MML_in$density_quant=factor(quant_conv[MML_in$density_quant],levels=quant_conv)
#PLotting boxplot
pdf(paste0(figure_path,'CpG_density_MML_features.pdf'),width=3.5,height=3.5)
ggplot(as.data.frame(mcols(MML_in)),aes(x=density_quant, y=MML))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()   

#Feature enrichment of low MML regions
genomic_features=readRDS(genomic_features_file)
#Figure S5
olap_islands=findOverlaps(MML_in,genomic_features$`CpG island`)
olap_shores=findOverlaps(MML_in,genomic_features$`CpG shore`)
olap_shelf=findOverlaps(MML_in,genomic_features$`CpG shelf`)
olap_open_sea=findOverlaps(MML_in,genomic_features$`CpG open sea`)

olap_islands=findOverlaps(MML_in,genomic_features$`CpG island`)

CpG_density_MML=rbind(data.table(MML=MML_in$MML[queryHits(olap_islands)],feature='islands'),
                      data.table(MML=MML_in$MML[queryHits(olap_shores)],feature='shores'),
                      data.table(MML=MML_in$MML[queryHits(olap_shelf)],feature='shelf'),
                      data.table(MML=MML_in$MML[queryHits(olap_open_sea)],feature='open sea'))

pdf(paste0(figure_path,'CpG_density_MML_features.pdf'),width=3.5,height=3.5)
ggplot(CpG_density_MML,aes(x=feature, y=MML))+
  ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
  ylab("MML")+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off() 

# In mouse context --------------------------------------------------------
CG_exp_agnostic_mm10_file='../downstream/output/mouse_analysis/CpG_density/analyzed_region_mm10.rds'
mml=readRDS(MML_matrix_file)
nme=readRDS(NME_matrix_file)
density_gr=unique(c(granges(nme),granges(mml)))
gr_seq=getSeq(Mmusculus,density_gr,as.character=T)
density_gr$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
saveRDS(density_gr,CG_exp_agnostic_mm10_file)
density_gr=readRDS(CG_exp_agnostic_mm10_file)
mm10_CpG=getCpgSitesmm10()

nme$CG_mm10=countOverlaps(nme,mm10_CpG)
mml$CG_mm10=countOverlaps(mml,mm10_CpG)
#the regions are the same for mml and nme so use same gr files

nme_olap=findOverlaps(nme,density_gr,type='equal')
nme$CGcont_exp[queryHits(nme_olap)]=density_gr$CGcont_exp[subjectHits(nme_olap)]
nme$density=nme$CG_mm10/nme$CGcont_exp
mml_olap=findOverlaps(mml,density_gr,type='equal')
mml$CGcont_exp[queryHits(mml_olap)]=density_gr$CGcont_exp[subjectHits(mml_olap)]
mml$density=mml$CG_mm10/mml$CGcont_exp
nme$density_quant=findInterval(nme$density,seq(0,1,0.1))
mml$density_quant=findInterval(mml$density,seq(0,1,0.1))
quant_conv=c(paste0(seq(0,0.9,0.1),'-',seq(0.1,1,0.1)),'>1')
nme$density_quant=factor(quant_conv[nme$density_quant],levels=quant_conv)

mml$density_quant=factor(quant_conv[mml$density_quant],levels=quant_conv)
density_mouse_calc<-function(gr_in,stat_name="NME"){
  gr_in=mcols(gr_in)
  gr_in=melt.data.table(as.data.table(gr_in),id.vars = c("CG_mm10","CGcont_exp","density","density_quant"),variable.name = "Sample",value.name="stat_in")
  gr_in$density_quant=findInterval(gr_in$density,seq(0,1,0.1))
  #NME_in$density_quant[NME_in$density_quant==6]=5#11th quantile is the maximum number, move to 10th
  quant_conv=c(paste0(seq(0,0.9,0.1),'-',seq(0.1,1,0.1)),'>1')
  gr_in$density_quant=factor(quant_conv[gr_in$density_quant],levels=quant_conv)
  print(ggplot(as.data.frame(gr_in),aes(x=density_quant, y=stat_in))+
          ylim(c(0,1))+geom_boxplot(outlier.shape = NA)+theme_glob+xlab("CpG density")+
          ylab(stat_name)+theme(axis.text.x =  element_text(angle = 90, vjust = 0.5, hjust=1)))
}
pdf(paste0(figure_path,'mouse_NME_density_boxplot.pdf'),width=3.5,height=3.5)
density_mouse_calc(nme,stat_name="NME")
dev.off()
pdf(paste0(figure_path,'mouse_MML_density_boxplot.pdf'),width=3.5,height=3.5)
density_mouse_calc(mml,stat_name="MML")
dev.off()

NME_dt=convert_GR(nme,direction='DT')
NME_dt_mt=melt.data.table(NME_dt,id.vars=c('region','density','CGcont_exp','CG_mm10'),value.name = "NME",variable.name='Sample')


cor.test(NME_dt_mt$NME,NME_dt_mt$density,method='pearson')

