rm(list=ls())
source("mainFunctions_sub.R")

theme_glob=theme_classic()+theme(plot.title = element_text(hjust = 0.5,size=10),
                                 axis.title.x=element_text(hjust=0.5,size=9,face="bold"),
                                 axis.title.y=element_text(hjust=0.5,size=9,face="bold"),
                                 axis.text.x=element_text(size=6),
                                 axis.text.y=element_text(size=6))
# get all variant to prepare motifbreakR analysis---------------------------------------------------------
variant_HetCpG=readRDS(variant_HetCpG_file)
names(variant_HetCpG)=NULL
saveRDS(unique(do.call('c',variant_HetCpG)),
        '../downstream/output/human_analysis/motif_analysis/variant_in_all.rds')

# Read in result ----------------------------------------------------------
#Checked in clean run 
motif_dir='../downstream/output/human_analysis/motif_analysis/JASPAR_default/'
motif_gene=GRanges()
for(fn in dir(motif_dir)){
  
  motif_gene=c(motif_gene,readRDS(paste0(motif_dir,fn)))
  
}
saveRDS(motif_gene,motif_gene_file)

# Find motif prefere NME or MML using binomial test using all the regions-----------------------
motif_gene <- readRDS(motif_gene_file)
variant_HetCpG_meta=readRDS(variant_HetCpG_meta_file)
  #NME

motif_dir_dNME=direction_calc_enriched_subj(motif_gene,variant_HetCpG_meta,
                                       unique(motif_gene$geneSymbol),pval_cutoff=0.1,stat="NME")
colnames(motif_dir_dNME)[c(1,6,7)]=c('TF','Pvalue','Proportion')
motif_dir_dNME$FDR=p.adjust(motif_dir_dNME$Pvalue,method="BH")
motif_dir_dNME$proportion_high_NME=motif_dir_dNME$Proportion
motif_dir_dNME$proportion_high_NME_CI=paste0("[",round(motif_dir_dNME$lowerCI,digits=2),' - ',round(motif_dir_dNME$upperCI,digits=2),"]")
motif_dir_dNME$proportion_low_NME=1-motif_dir_dNME$Proportion
motif_dir_dNME$proportion_low_NME_CI=paste0("[",round(1-motif_dir_dNME$upperCI,digits=2),' - ',round(1-motif_dir_dNME$lowerCI,digits=2),"]")

saveRDS(motif_dir_dNME,'../downstream/output/human_analysis/motif_analysis/dNME_all.rds')
motif_dir_dNME=readRDS('../downstream/output/human_analysis/motif_analysis/dNME_all.rds')
write.csv(motif_dir_dNME[FDR<=0.1&Proportion>0.5,list(TF,proportion_high_NME,Pvalue,FDR,proportion_high_NME_CI)], row.names =F,
          '../downstream/output/graphs_tables/motif_preference_table/All_regions/table1_motif_prefer_high_NME.csv')
write.csv(motif_dir_dNME[FDR<=0.1&Proportion<0.5,list(TF,Proportion=proportion_low_NME,Pvalue,FDR,proportion_low_NME_CI)], row.names =F,
          '../downstream/output/graphs_tables/motif_preference_table/All_regions/table2_motif_prefer_low_NME.csv')
write.table(unlist(strsplit(gsub('\\(.*',"",motif_dir_dNME[order(FDR,-Proportion )]$TF),"::")),'../downstream/output/human_analysis/motif_analysis/dNME_bg_gene.txt',row.names = F,col.names = F,quote=F)
NME_only_GO=GO_run(unlist(strsplit(gsub('\\(.*',"",motif_dir_dNME[FDR<=0.1]$TF),"::")),
                   unlist(strsplit(gsub('\\(.*',"",motif_dir_dNME[order(FDR,-Proportion )]$TF),"::")),1,mapping="org.Hs.eg.db")
write.csv(NME_only_GO[FDR<=0.2&FC>=1.5],'../downstream/output/human_anacelysis/motif_analysis/NME_only_topGO.csv')
pdf('../downstream/output/human_analysis/motif_analysis/ASCL1_NME_allele.pdf',width=2.35,height=2.35)
motif_out_ASCL1=plot_merge_SNP_motif(variant_HetCpG_meta,motif_gene,motif="ASCL1",stat="NME",pval_cutoff=pval_cutoff,theme_glob=theme_glob)
dev.off()
#MML
motif_dir_dMML=direction_calc_enriched_subj(motif_gene,variant_HetCpG_meta,
                                       unique(motif_gene$geneSymbol),pval_cutoff=0.1,stat="MML")
colnames(motif_dir_dMML)[c(1,6,7)]=c('TF','Pvalue','Proportion')
motif_dir_dMML$FDR=p.adjust(motif_dir_dMML$Pvalue,method="BH")
motif_dir_dMML$Proportion_high_MML=motif_dir_dMML$Proportion
motif_dir_dMML$proportion_high_MML_CI=paste0("[",round(motif_dir_dMML$lowerCI,digits=2),' - ',round(motif_dir_dMML$upperCI,digits=2),"]")
motif_dir_dMML$Proportion_low_MML=1-motif_dir_dMML$Proportion_high_MML
motif_dir_dMML$proportion_low_MML_CI=paste0("[",round(1-motif_dir_dMML$upperCI,digits=2),' - ',round(1-motif_dir_dMML$lowerCI,digits=2),"]")
saveRDS(motif_dir_dMML,'../downstream/output/human_analysis/motif_analysis/dMML_all.rds')
motif_dir_dMML=readRDS('../downstream/output/human_analysis/motif_analysis/dMML_all.rds')
write.csv(motif_dir_dMML[FDR<=0.1&Proportion<0.5,list(TF,Proportion=Proportion_low_MML,Pvalue,FDR,proportion_low_MML_CI)], row.names =F,
          '../downstream/output/graphs_tables/motif_preference_table/All_regions/motif_prefer_low_MML.csv')
write.csv(motif_dir_dMML[FDR<=0.1&Proportion>0.5,list(TF,Proportion_high_MML,Pvalue,FDR,proportion_high_MML_CI)], row.names =F,
          '../downstream/output/graphs_tables/motif_preference_table/All_regions/motif_prefer_high_MML.csv')
write.table(unlist(strsplit(gsub('\\(.*',"",motif_dir_dMML[order(FDR,-Proportion )]$TF),"::")),'../downstream/output/human_analysis/motif_analysis/dMML_bg_gene.txt',row.names = F,col.names = F,quote=F)
MML_only_GO=GO_run(unlist(strsplit(gsub('\\(.*',"",motif_dir_dMML[FDR<=0.1]$TF),"::")),
                   unlist(strsplit(gsub('\\(.*',"",motif_dir_dMML[order(FDR,-Proportion )]$TF),"::")),1,mapping="org.Hs.eg.db")
write.csv(MML_only_GO[FDR<=0.2&FC>=1.5],'../downstream/output/human_analysis/motif_analysis/MML_only_topGO.csv')

pdf('../downstream/output/human_analysis/motif_analysis/LEF1_MML_allele.pdf',width=2.35,height=2.35)
plot_merge_SNP_motif(variant_HetCpG_meta,motif_gene,motif="LEF1",stat="MML",pval_cutoff=pval_cutoff,theme_glob=theme_glob)
dev.off()
# Find low MML and high NME only ------------------------------------------
#Rename columns
low_MML=motif_dir_dMML[FDR<=0.1&Proportion<0.5]$TF
high_NME=motif_dir_dNME[FDR<=0.1&Proportion>0.5]$TF
motif_dir_dMML$dMML_pvalue=motif_dir_dMML$Pvalue
motif_dir_dMML$dMML_FDR=motif_dir_dMML$FDR
motif_dir_dMML$CI=motif_dir_dMML$proportion_low_MML_CI
motif_dir_dNME$dNME_pvalue=motif_dir_dNME$Pvalue
motif_dir_dNME$dNME_FDR=motif_dir_dNME$FDR
motif_dir_dNME$CI=motif_dir_dNME$proportion_high_NME_CI
#Find low MML only and high NME only
low_MML_TF=low_MML[!low_MML%in%high_NME]
low_MML_only=cbind(motif_dir_dMML[match(low_MML_TF,TF)],
                   motif_dir_dNME[match(low_MML_TF,TF)])
high_NME_TF=high_NME[!high_NME%in%low_MML]
high_NME_only=cbind(motif_dir_dNME[match(high_NME_TF,TF)],
                    motif_dir_dMML[match(high_NME_TF,TF)])

write.csv(low_MML_only[order(Proportion_low_MML),list(TF,Proportion_low_MML,dMML_pvalue,dMML_FDR,proportion_high_NME,dNME_pvalue,dNME_FDR,CI)],
          '../downstream/output/graphs_tables/motif_preference_table/All_regions/motif_prefer_low_MML_only.csv',row.names = F)
write.csv(high_NME_only[order(proportion_high_NME),list(TF,proportion_high_NME,dNME_pvalue,dNME_FDR,Proportion_low_MML,dMML_pvalue,dMML_FDR,CI)],
          '../downstream/output/graphs_tables/motif_preference_table/All_regions/motif_prefer_high_NME_only.csv',row.names = F)

# find OMIM annotations ---------------------------------------------------
low_MML_motif=fread('../downstream/output/graphs_tables/motif_preference_table/All_regions/motif_prefer_low_MML_only.csv')
high_NME_motif=fread('../downstream/output/graphs_tables/motif_preference_table/All_regions//motif_prefer_high_NME_only.csv')
#Requested from:https://www.omim.org/downloads
OMIM=fread('../downstream/input/human_analysis/SNP_biology/genemap2.txt',skip=3)
OMIM$`Gene Symbols`=strsplit(as.character(OMIM$`Gene Symbols`),', ')
OMIM=OMIM[Phenotypes!=""]
motif_ent_OMIM=OMIM_annotation(high_NME_motif,OMIM)
motif_ent_only_OMIM=OMIM_annotation(high_NME_motif[!(TF%in%low_MML_motif$TF)],OMIM)
motif_low_MML_only_OMIM=OMIM_annotation(low_MML_motif[!(TF%in%high_NME_motif$TF)],OMIM)
motif_shared_OMIM=OMIM_annotation(low_MML_motif[TF%in%high_NME_motif$TF],OMIM)
write.csv(motif_low_MML_only_OMIM, row.names =F,
          "../downstream/output/human_analysis/motif_analysis/motif_prefer_low_MML_only_OMIM.csv")
write.csv(motif_ent_only_OMIM, row.names =F,
          "../downstream/output/human_analysis/motif_analysis/motif_prefer_high_NME_only_OMIM.csv")


