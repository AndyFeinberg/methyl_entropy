rm(list=ls())
source("mainFunctions_sub.R")
# Human analysis ----------------------------------------------------------
#ENCFF001TDO from ENCODE is the same file
#https://github.com/Boyle-Lab/Blacklist/blob/master/lists/Blacklist_v1/
hg19_bl=import.bed('../downstream/input/human_analysis/QC/hg19_blacklist.bed.gz')
#Get all CpG location
chrs <- names(Hsapiens)[1:22]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:22, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 1)))) #use first location
# 20kb from hg19 TSS ------------------------------------------------------
extend_size=20000
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- GenomicFeatures::genes(txdb)
TSS<-promoters(genes,upstream=100,downstream=100)
TSS$gene_symbol=AnnotationDbi::select(Homo.sapiens,key=as.character(TSS$gene_id),
                                      keytype="ENTREZID",columns=c("SYMBOL"))$SYMBOL
TSS=TSS[seqnames(TSS) %in% paste0('chr',1:22)]
start(TSS)=start(TSS)-extend_size
end(TSS)=end(TSS)+extend_size+1
TSS=trim(TSS)
TSS_break=subdivideGRanges(TSS,250)
gff_gen(TSS_break,cpgr,hg19_bl,'../downstream/output/human_analysis/CPEL_inputs/hg19_TSS_20kb_250bp.rds',
        '../downstream/output/human_analysis/CPEL_inputs/hg19_TSS_20kb_250bp.gff')

# DNase vs non-DNase ------------------------------------------------------
#Those are from ken
DNase=readRDS('../downstream/input/human_analysis/DNase_hg19_250bp.rds')
control=readRDS('../downstream/input/human_analysis/DNase_hg19_250bp_control.rds')
TSS_break=c(DNase,control)
gff_gen(TSS_break,cpgr,hg19_bl,'../downstream/output/human_analysis/CPEL_inputs/hg19_DNase_250bp_gff.rds','../downstream/output/human_analysis/CPEL_inputs/hg19_DNase_250bp_gff.gff')

# All region have SNP -----------------------------------------------------
GR_merge=readRDS(GR_merge_file)
gff_gen(granges(unique(GR_merge)),cpgr,hg19_bl,'../downstream/output/human_analysis/CPEL_inputs/human_ASM_region_allele_agnostic_250bp.rds',
        '../downstream/output/human_analysis/CPEL_inputs/human_ASM_region_allele_agnostic_250bp.gff')


# rest of the regions in hg19 that has not been analyzed ------------------
hg19_all=GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg19))
hg19_all=hg19_all[seqnames(hg19_all) %in% paste0('chr',c(1:22,'X','Y'))]
analyzed_regions=c(readRDS('../downstream/input//human_analysis/DNase_hg19_250bp.rds'),
                   readRDS('../downstream/input/human_analysis/DNase_hg19_250bp_control.rds'),
                   readRDS('../downstream/output/human_analysis/CPEL_inputs/hg19_TSS_20kb_250bp.rds'))
analyzed_regions=c(granges(analyzed_regions),granges(unique(readRDS(GR_merge_file))))
hg19_all_comp=setdiff(hg19_all,analyzed_regions)
hg19_all_comp_break=subdivideGRanges(hg19_all_comp,250)
chrs <- names(Hsapiens)[1:24]#2276796
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 1)))) #use first location


gff_gen(hg19_all_comp_break,cpgr,hg19_bl,'../downstream/output/human_analysis/CPEL_inputs/hg19_allele_agnostic_analysis_compliment.rds',
        '../downstream/output/human_analysis/CPEL_inputs/hg19_allele_agnostic_analysis_comp.gff')

# Linux code for correcting format and copy to each file ------------------
#Running dir: /ibox/afeinbe2/yfang/allele_specific_roadmap_CEPL/work_archive/CpelAsm/
# sed -i 's/%2c/,/g' hg19_TSS_20kb_250bp.gff
# sed -i 's/chr//g' hg19_TSS_20kb_250bp.gff
# subj=(H9 HUES64 skin03 STL001 STL002 STL003 STL011 H1 HuFGM02 112 149 150)
# for i in "${subj[@]}"; do cp hg19_TSS_20kb_250bp.gff ~/work/shared/CpelAsm/data/${i}/cpelasm/${i}_allele_agnostic_analysis.gff; done
# sed -i 's/%2c/,/g' hg19_allele_agnostic_analysis_comp.gff
# sed -i 's/];/]/g' hg19_allele_agnostic_analysis_comp.gff
# sed -i 's/chr//g' hg19_allele_agnostic_analysis_comp.gff
# subj=(H9 HUES64 skin03 STL001 STL002 STL003 STL011 H1 HuFGM02 112 149 150)
# for i in "${subj[@]}"; do cp hg19_allele_agnostic_analysis_comp.gff /ibox/afeinbe2/yfang/allele_specific_roadmap_CEPL/work_archive/CpelAsm/data/${i}/cpelasm/${i}_allele_agnostic_analysis.gff; done
#Clean human analysis
rm(list=ls())
source("mainFunctions_sub.R")
# mm10 allele-agnostic analysis for DNase-control, PRC2 binding regions -------------------------------------------
#https://github.com/Boyle-Lab/Blacklist/blob/master/lists/Blacklist_v1/
blacklist_region=import.bed('../downstream/input/mouse_analysis/mm10.blacklist.bed.gz')
chrs <- names(Mmusculus)[1:21]#2276796
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Mmusculus[[x]])))
cpgr <- do.call(c, lapply(1:21, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs[[x]], width = 1)))) #use first location
#DNase vs control is also from Ken
DNAase=readRDS('../downstream/input/mouse_analysis/DNase_mm10_peak_merge_250bp.rds')
control=readRDS('../downstream/input/mouse_analysis/DNase_mm10_peak_merge_250bp_control.rds')
#PRC is from this command: takes long time, from Bin's paper
#pooled_PRC= read_chromHMM_bed_PRC('../downstream/input/mouse_analysis/enhancer_selection/chromHMM_mm10/pooled/','pooled',blacklist_region,cpgr)
PRC_binding=readRDS('../downstream/input/mouse_analysis/mm10_PRC.rds')
PRC_binding=granges(PRC_binding)


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- GenomicFeatures::genes(txdb)
TSS<-promoters(genes,upstream=0,downstream=0)
TSS$gene_name=AnnotationDbi::select(Mus.musculus,key=as.character(TSS$gene_id),
                                    keytype="ENTREZID",columns=c("SYMBOL"))$SYMBOL

PRC_binding$region_type="PRC"
DNAase$region_type='DNase'
control$region_type='control'
TSS_break=c(DNAase,control,PRC_binding)
#dis com
TSS_break=dist_calc(TSS_break,TSS)
gff_gen(TSS_break,cpgr,blacklist_region,'../downstream/output/mouse_analysis/CPEL_inputs/mm10_allele_agnostic_analysis.rds',
        mouse_DNase_control_gff_file)
# mm10 allele-agnostic analysis for complementary regions -------------------------------------------
mm10_all=GRanges(seqinfo(BSgenome.Mmusculus.UCSC.mm10))
mm10_all=mm10_all[seqnames(mm10_all) %in% paste0('chr',c(1:21,'X','Y'))]
analyzed_regions=readGFFAsGRanges(mouse_DNase_control_gff_file)
mm10_all_comp=setdiff(mm10_all,analyzed_regions)
mm10_all_comp_break=subdivideGRanges(mm10_all_comp,250)
chrs <- names(Mmusculus)[1:21]#2276796
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Mmusculus[[x]])))
cpgr <- do.call(c, lapply(1:21, function(x) GRanges(names(Mmusculus)[x], IRanges(cgs[[x]], width = 1)))) #use first location
gff_gen(mm10_all_comp_break,cpgr,blacklist_region,'../downstream/output/mouse_analysis/mm10_allele_agnostic_analysis_compliment.rds',
        mouse_compliment_gff_file)
#Sanity check
seqlength=0
for(i in c(1:19,'X','Y')){seqlength=seqlength+length(BSgenome.Mmusculus.UCSC.mm10[[paste0('chr',i)]])} 
mm10_all_break=subdivideGRanges(mm10_all,250)
mm10_all_break$CpG=countOverlaps(mm10_all_break,cpgr)
sum(width(mm10_all_break[mm10_all_break$CpG>0]))/seqlength#0.7355
# mm10 linux reformat -----------------------------------------------------
#mm10
# sed -i 's/%2c/,/g' mm10_allele_agnostic_analysis.gff
# sed -i 's/];$/]/g' mm10_allele_agnostic_analysis.gff
# sed -i 's/rtracklayer/\./g' mm10_allele_agnostic_analysis.gff
# sed -i 's/sequence_feature/\./g' mm10_allele_agnostic_analysis.gff


  