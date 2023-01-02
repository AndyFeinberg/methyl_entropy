#!/usr/bin/env Rscript
#Nedd R/3.6.1, gcc/5.5.0
source('motifbreakR_parallel.R')
setwd("../")
args = commandArgs(trailingOnly=TRUE)
##############################Motifbreak_R analysis###################################
if (!requireNamespace("MotifDb", quietly = TRUE))
{BiocManager::install("MotifDb")}
library(MotifDb)
#need libatlas for linux: ml atlas
if (!requireNamespace("motifbreakR", quietly = TRUE))
{BiocManager::install("motifbreakR")}
library(motifbreakR)
if (!requireNamespace("SNPlocs.Hsapiens.dbSNP142.GRCh37", quietly = TRUE))
{BiocManager::install("SNPlocs.Hsapiens.dbSNP142.GRCh37")}
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE))
{BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")}
library("BSgenome.Hsapiens.UCSC.hg19")
if (!requireNamespace("BiocParallel", quietly = TRUE))
{BiocManager::install("BiocParallel")}
library(BiocParallel)
library('org.Hs.eg.db')
library(VGAM)
if (!requireNamespace("TFBSTools", quietly = TRUE))
{ devtools::install_github("ge11232002/TFBSTools")}
library(TFBSTools)
motif_break<-function(gr_in,motif_list){
  gr_in_gr=granges(gr_in)
  strand(gr_in_gr)='*'
  attributes(gr_in_gr)$genome.package="BSgenome.Hsapiens.UCSC.hg19"
  #Make sure ref and alt agree with BS genome #may not necessary
  ref_BS=as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,gr_in_gr))
  alt_BS=gr_in$ALT
  switch_nucleotide=which(gr_in$REF!=ref_BS)
  print(switch_nucleotide)
  alt_BS[switch_nucleotide]=gr_in$REF[switch_nucleotide]
  #Get the necessary information for program running
  gr_in_gr$REF=unlist(DNAStringSetList(ref_BS))
  gr_in_gr$ALT=unlist(DNAStringSetList(alt_BS))
  names(gr_in_gr)=gr_in$snpId
  print(gr_in_gr)

  return(gr_in_gr)
}
n=as.numeric(args[1])
file_in=args[2]
cut_number=as.numeric(args[3])
file_out=args[4]
method_motif=args[5]
cat('Reading in',file_in,'for: ',n,'with gap of',1/cut_number,'\n')
tt1=proc.time()[[3]]
variant_in=readRDS(paste('../downstream/output/human_analysis/motif_analysis/',file_in,sep=''))
JASPAR_2020=readRDS('../downstream/input/human_analysis/motif_analysis/JASPAR_2020_human_PFM.rds')
meta_JASPAR=data.frame(providerName=names(JASPAR_2020),providerId=unlist(lapply(JASPAR_2020@listData,function(x) x@ID)),
                       datasource='JASPAR2020',geneSymbol=unlist(lapply(JASPAR_2020@listData,function(x) x@name)),stringsAsFactors = F)
meta_JASPAR$geneId=mapIds(org.Hs.eg.db,meta_JASPAR$geneSymbol,'ENTREZID','SYMBOL')
meta_JASPAR$geneIdType=NA
meta_JASPAR$geneIdType[!is.na(meta_JASPAR$geneId)]='ENTREZID'
meta_JASPAR$organism='Hsapiens'
meta_JASPAR[,c('protein','proteinIdType','sequenceCount','bingdingSequence','bindingDomain','tfFamily','experiment','pubmedID')] =NA
mcols(JASPAR_2020)=meta_JASPAR
JASPAR_2020@listData=lapply(JASPAR_2020@listData,function(x) x@profileMatrix)
#function to separate them into  non-overlapping ?: split
split_data=cut(1:length(variant_in),breaks=cut_number,label=FALSE)
tt2=proc.time()[[3]]
cat('Finishing reading in data in: ',tt2-tt1,'\n')
cat('Starting motifbreak\n')
motif_out=motif_break(variant_in[split_data==n],JASPAR_2020)
cat('Starting motifbreakR\n')
#HS.SELEX=subset (MotifDb, organism=='Hsapiens'&dataSource=="jolma2013")
#use 2 statistic: try ic
motif_out<- motifbreakR_parallel(snpList = motif_out, filterp = TRUE,
                       pwmList = JASPAR_2020,
                       threshold = 1e-4,
                       method = method_motif,
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),verbose=T,
                       BPPARAM =MulticoreParam(workers=22,progressbar = TRUE))

tt3=proc.time()[[3]]
cat('Motif break finish in',tt3-tt2,'\n')
saveRDS(motif_out,paste('../downstream/output/human_analysis/motif_analysis/JASPAR_default/',file_out,'_',n,'_',method_motif,'.rds',sep=''))
tt4=proc.time()[[3]]
cat('saving data finish in',tt4-tt3,'\n')



