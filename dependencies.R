suppressPackageStartupMessages({
  if (!requireNamespace("ggfortify", quietly = TRUE))
  {
    install.packages("glue")
    BiocManager::install("ggfortify")
    }
  library(ggfortify)
  if (!requireNamespace("xlsx", quietly = TRUE))
  {
    install.packages("xlsx")
  }
  library(xlsx)
  if (!requireNamespace("shiny", quietly = TRUE))
  {install.packages("htmltools")
    BiocManager::install("shiny")}
  library(shiny)
  if (!requireNamespace("RColorBrewer", quietly = TRUE))
  {BiocManager::install("RColorBrewer")}
  library(RColorBrewer)
  if (!requireNamespace("ggplot2", quietly = TRUE))
  {BiocManager::install("ggplot2")}
  library(ggplot2)
  # if (!requireNamespace("ggpubr", quietly = TRUE))
  # {BiocManager::install("ggpubr")}
  # library(ggpubr)
  
  if (!requireNamespace("lattice", quietly = TRUE))
  {install.packages("lattice")}
  library(lattice)
  if (!requireNamespace("gridExtra", quietly = TRUE))
  {BiocManager::install("gridExtra")}
  library(gridExtra)
  if (!requireNamespace("grid", quietly = TRUE))
  {BiocManager::install("grid")}
  library(grid)
  if (!requireNamespace("pander", quietly = TRUE))
  {install.packages("pander")}
  library(pander)
 if (!requireNamespace("preprocessCore", quietly = TRUE))
    BiocManager::install("preprocessCore")
  library(preprocessCore)
  if (!requireNamespace("rtracklayer", quietly = TRUE))
  {BiocManager::install("rtracklayer")}
  library(rtracklayer)
  if (!requireNamespace("GenomicRanges", quietly = TRUE))
  {BiocManager::install("GenomicRanges")}
  library(GenomicRanges)
  if (!requireNamespace("data.table", quietly = TRUE))
  {BiocManager::install("data.table")}
  library(data.table)
  if (!requireNamespace("Gmisc", quietly = TRUE))
  {BiocManager::install("Gmisc")}
  library(Gmisc)
  if (!requireNamespace("matrixStats", quietly = TRUE))
  {BiocManager::install("matrixStats")}
  library(matrixStats)
  if (!requireNamespace("parallel", quietly = TRUE))
  {BiocManager::install("parallel")}
  library(parallel)
  if (!requireNamespace("readr", quietly = TRUE))
  {install.packages("readr")}
  library(readr)
  if (!requireNamespace("readxl", quietly = TRUE))
  {install.packages("readxl")}
  library(readxl)
  if (!requireNamespace("reshape2", quietly = TRUE))
  {BiocManager::install("reshape2")}
  library(reshape2)
  if (!requireNamespace("exomeCopy", quietly = TRUE))
  {BiocManager::install("exomeCopy")}
  library(exomeCopy)
  if (!requireNamespace("pheatmap", quietly = TRUE))
  {BiocManager::install("pheatmap")}
  library(pheatmap)
  #hg19 genome 
  if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE))
  {BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")}
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  if (!requireNamespace("Homo.sapiens", quietly = TRUE))
  {BiocManager::install("Homo.sapiens")}
  library(Homo.sapiens)
  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE))
  {BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")}
  library(BSgenome.Hsapiens.UCSC.hg19)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  {BiocManager::install("org.Hs.eg.db")}
  library(org.Hs.eg.db)
  
  #hg38 genome: not used maybe
  if (!requireNamespace("SNPlocs.Hsapiens.dbSNP151.GRCh38", quietly = TRUE)){
    BiocManager::install("SNPlocs.Hsapiens.dbSNP151.GRCh38")}
  library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
  
  #mm10 genome related
  if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)){
    BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")}
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  if (!requireNamespace("Mus.musculus", quietly = TRUE)){
    BiocManager::install("Mus.musculus")}
  library(Mus.musculus)
  if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)){
    BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")}
  library(BSgenome.Mmusculus.UCSC.mm10)
  
  #genome tools
  if (!requireNamespace("liftOver", quietly = TRUE)){
    BiocManager::install("liftOver")}
  library(liftOver)
  
  #Accessibility analysis
  if (!requireNamespace("limma", quietly = TRUE))
  {BiocManager::install("limma")}
  library(limma)
  if (!requireNamespace("edgeR", quietly = TRUE))
  {BiocManager::install("edgeR")}
  library(edgeR)
  if (!requireNamespace("tximport", quietly = TRUE))
  {BiocManager::install("tximport")}
  library(tximport)
  if (!requireNamespace("DESeq2", quietly = TRUE))
  {BiocManager::install("DESeq2")}
  library(DESeq2)
  
  #Annotations
  if (!requireNamespace("annotatr", quietly = TRUE))
  {BiocManager::install("annotatr",INSTALL_opts = c('--no-lock'))}
  library(annotatr)
  if (!requireNamespace("AnnotationHub", quietly = TRUE))
  {BiocManager::install("AnnotationHub")}
  library(AnnotationHub)
  if (!requireNamespace("GenomicFeatures", quietly = TRUE))
  {BiocManager::install("GenomicFeatures")}
  library(GenomicFeatures)
  if (!requireNamespace("VariantAnnotation", quietly = TRUE))
  {BiocManager::install("VariantAnnotation")}
  library(VariantAnnotation)
  if (!requireNamespace("topGO", quietly = TRUE))
  {BiocManager::install("topGO")}
  library(topGO)
  if (!requireNamespace("ggpubr", quietly = TRUE))
  {install.packages("ggpubr")}
  library(ggpubr)
  #Motif analysis: need atlas, gcc/5.5.0, R/3.6.1
  if (!requireNamespace("motifbreakR", quietly = TRUE)){
    BiocManager::install("motifbreakR")}
  library(motifbreakR)
  if (!requireNamespace("TFBSTools", quietly = TRUE)){
    BiocManager::install("TFBSTools")}
  library(TFBSTools)
  if (!requireNamespace("MotifDb", quietly = TRUE))
  {BiocManager::install("MotifDb")}
  library(MotifDb)
  #correlation
  
  if (!requireNamespace("MASS", quietly = TRUE)){
   install.packages("MASS")}
  library(MASS)
  


#Pubmed
  if (!requireNamespace("rentrez", quietly = TRUE)){
   install.packages("rentrez")}
  library(rentrez)
  if (!requireNamespace("easyPubMed", quietly = TRUE)){
   install.packages("easyPubMed")}
  library(easyPubMed)
})
# Not in use --------------------------------------------------------------
# 
# if (!requireNamespace("Repitools", quietly = TRUE))
# {BiocManager::install("Repitools")}
# library(Repitools)
# if(!require(psych)){install.packages("psych")}
# library('psych')
# if(!require(vcd)){install.packages("vcd")}
# if(!require(DescTools)){install.packages("DescTools")}
# if(!require(rcompanion)){install.packages("rcompanion")}
# if(!require(ggpubr)){install.packages("ggpubr")}
# library("ggpubr")
# if(!require(readxl)){install.packages("readxl")}
# library(readxl) 
# if(!require(ggpointdensity)){install.packages("ggpointdensity")}
# library(ggpointdensity)
# if (!requireNamespace("rethinking", quietly = TRUE))
# {install.packages(c("coda","mvtnorm","devtools","loo","dagitty"))
#   library(devtools)
#   devtools::install_github("rmcelreath/rethinking")}
# library(rethinking)
# if (!requireNamespace("gwasrapidd", quietly = TRUE)){
#   remotes::install_github("ramiromagno/gwasrapidd")}
# library(gwasrapidd)
# library(qusage)
# if (!requireNamespace("gwascat", quietly = TRUE)){
# BiocManager::install("gwascat")}
# library(gwascat)
