# Dependencies in use -----------------------------------------------------
#ml java/13.0.2
#Need R/3.6.1, gcc/5.5.0, and atlas
#clean run: remove: rm -r /home-net/home-4/yfang27@jhu.edu/R/x86_64-pc-linux-gnu-library/3.6/gcc/5.5/*
#rm -r /scratch/users/yfang27@jhu.edu/yfang/temp_all/*
#plot 
# get all file location ---------------------------------------------------
source('file_path.R')
source('dependencies.R')
#Set DT use only 1 thread
setDTthreads(1)
# main functions ----------------------------------------------------------
#Get CpG sites from hg19
getCpgSitesH19 <- function(chrsOfInterest=paste("chr",c(1:22,"X","Y"),sep="")){
  # Obtain all CpG sites
  cgs <- lapply(chrsOfInterest, function(x) GRanges(x,IRanges(start(matchPattern("CG", Hsapiens[[x]])),width=2)))
  
  # Set genome and seqlengths
  cgs <- setGenomeLengths(do.call('c',cgs))
  # Return CpG site GR
  return(cgs)
}
#Color theme for mouse
mouse_color<-function(){
  tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
  cv <- brewer.pal(length(tissue_all),'Set1')
  cv[6] <- 'goldenrod1'
  names(cv) <-tissue_all
  return(cv)
}
#generating gff file
gff_gen<-function(TSS_break,cpgr,blacklist_region,rds_save_file,out_name){
  olap=findOverlaps(TSS_break,cpgr)
  CpG_df=data.frame(TSS_hit=queryHits(olap),CpG_start=start(cpgr)[subjectHits(olap)])
  CpG_df_agg=aggregate(CpG_start~TSS_hit,CpG_df,function(x) paste(x,collapse=', '))
  CpG_df_agg$N=unlist(lapply(CpG_df_agg$CpG_start,function(x) length(strsplit(x,', ')[[1]])))
  CpG_df_agg$CpG_start=paste("[",CpG_df_agg[,2],"]",sep='')
  TSS_break_out=TSS_break[CpG_df_agg[,1]]
  strand(TSS_break_out)='*'
  TSS_break_out$N=CpG_df_agg$N
  TSS_break_out$CpGs=CpG_df_agg$CpG_start
  #filter out blacklist region and region with N=1
  
  #TSS_break_out=TSS_break_out[TSS_break_out$N>1]
  olap=findOverlaps(TSS_break_out,blacklist_region)
  TSS_break_out=TSS_break_out[-queryHits(olap)]
  TSS_break_out_gff=granges(TSS_break_out)
  mcols(TSS_break_out_gff)=mcols(TSS_break_out)[,c('N','CpGs')]

  export.gff3(sort(unique(TSS_break_out_gff)),out_name)
  saveRDS(TSS_break_out,rds_save_file)
}
#Checking chromosome name
chr_check<-function(gr_in){
  #Get plus one location
  chr_check=grepl('chr',seqlevels(gr_in))#checking if seqlevels contain chr
  if(any(!chr_check)){seqlevels(gr_in)[which(!chr_check)]=paste("chr",seqlevels(gr_in)[which(!chr_check)],sep='')}
  return(gr_in)
}
#Read in gff file
readAllGff <- function(inDir,subjects,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Loop over all GFF files
  outGR <- GRanges()
  #subjects <- c("H9","HUES64","skin03","HuFGM02","STL001","STL002","STL003","STL011")
  for (sub in subjects) {
    # Import GFF file
    cat('importing',sub,'\n')
    tmpGR <- readSubGff(inDir,sub,chrsOfInterest)
    start(tmpGR)=start(tmpGR)
    # Append to output GR
    outGR <- c(outGR,tmpGR)
  }
  
  # Return GR with all haplotypes
  return(outGR)
  
}
#Read in gff file for each subject
readSubGff <- function(inDir,sub,chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Import GFF
  outGR <- import.gff(paste(inDir,sub,"_het.cpelasm.gff",sep=""))
  
  # Retain a the required subset
  outGR <- outGR[,c("N","CpGs","hetCpGg1","hetCpGg2")]
  
  # Add subject metadata column
  outGR$Subject <- sub
  
  # Make N numeric and filter
  outGR$N <- as.numeric(outGR$N)
  
  # Add genome info 
  outGR <- setGenomeLengths(outGR,chrsOfInterest=chrsOfInterest)
  
  # Return GR with all haplotypes
  return(outGR)
  
}
###Here most current version of Het CpG analysis

#From vcf file, extract het CpG information
extractHetCpG<-function(vcfDir,sub,chrsOfInterest=paste("chr",1:22,sep="")){
  cat('Processing subject:', sub,'\n')
  tt1=proc.time()[[3]]
  vcf <- readVcf(file=paste(vcfDir,sub,".phased.vcf.gz",sep=""),"hg19")
  
  gt <- as.vector(geno(vcf)$GT)
  vcf <- rowRanges(vcf)
  vcf$GT <- gt
  vcf$snpId <- paste(sub,seq(1:length(vcf)),sep="-")
  # Keep only relevant variables
  vcf <- vcf[,c("REF","ALT","GT","snpId")]
  vcf$REF <- as.character(vcf$REF)
  vcf$ALT <- as.character(unlist(vcf$ALT))
  names(vcf)=NULL
  # Delete labels
  vcf=het_CpG_df(vcf)
  vcf$sub=sub
  cat('Done processing',sub,'using:', proc.time()[[3]]-tt1,'\n')
  return(vcf)
  
}
#ID heterozygous CpG from vcf file
het_CpG_df<-function(var_gr){
  var_gr=chr_check(var_gr)
  plus_loc=as.character(Views(Hsapiens,GenomicRanges::shift(var_gr,1)))
  minus_loc=as.character(Views(Hsapiens,GenomicRanges::shift(var_gr,-1)))
  #get dinucleotide for ref, alt, plus and minus, find some examples region randomly: check if those match
  var_gr$REF_plus=paste0(var_gr$REF,plus_loc)
  var_gr$REF_minus=paste0(minus_loc,var_gr$REF)
  var_gr$ALT_plus=paste0(var_gr$ALT,plus_loc)
  var_gr$ALT_minus=paste0(minus_loc,var_gr$ALT)
  #get trinucleotide
  var_gr$REF_tri=paste0(minus_loc,var_gr$REF,plus_loc)
  var_gr$ALT_tri=paste0(minus_loc,var_gr$ALT,plus_loc)
  #currently not in use
  # #check if heterogygouze: note rowSum =2 have trinucleotide form CGG with ref =G alt =C
  # var_gr$npmCG=rowSums(as.data.table(mcols(var_gr))
  #                      [,.(str_count(REF_tri,pattern="CG"),str_count(ALT_tri,pattern="CG"))])
  # var_gr$HetCpG=var_gr$npmCG>0
  #Add CG content for genome1 and genome2 based on GT
  var_gr$genome1_plus=NA
  var_gr$genome1_minus=NA
  var_gr$genome1_tri=NA
  var_gr$genome2_plus=NA
  var_gr$genome2_minus=NA
  var_gr$genome2_tri=NA
  ####Fix the issue with genome file order does not equal to the ref/alt order, calculate the CG or het CG in genome 1 based on actual ref/alt order
  #Genome1
  var_gr$genome1_plus[var_gr$GT %in% c("0/1","0|1")]=var_gr$REF_plus[var_gr$GT %in% c("0/1","0|1")]
  var_gr$genome1_minus[var_gr$GT %in% c("0/1","0|1")]=var_gr$REF_minus[var_gr$GT %in% c("0/1","0|1")]
  var_gr$genome1_tri[var_gr$GT %in% c("0/1","0|1")]=var_gr$REF_tri[var_gr$GT %in% c("0/1","0|1")]
  
  var_gr$genome1_plus[var_gr$GT %in% c("1/0","1|0")]=var_gr$ALT_plus[var_gr$GT %in% c("1/0","1|0")]
  var_gr$genome1_minus[var_gr$GT %in% c("1/0","1|0")]=var_gr$ALT_minus[var_gr$GT %in% c("1/0","1|0")]
  var_gr$genome1_tri[var_gr$GT %in% c("1/0","1|0")]=var_gr$ALT_tri[var_gr$GT %in% c("1/0","1|0")]
  #Genome2
  var_gr$genome2_plus[var_gr$GT %in% c("0/1","0|1")]=var_gr$ALT_plus[var_gr$GT %in% c("0/1","0|1")]
  var_gr$genome2_minus[var_gr$GT %in% c("0/1","0|1")]=var_gr$ALT_minus[var_gr$GT %in% c("0/1","0|1")]
  var_gr$genome2_tri[var_gr$GT %in% c("0/1","0|1")]=var_gr$ALT_tri[var_gr$GT %in% c("0/1","0|1")]
  
  var_gr$genome2_plus[var_gr$GT %in% c("1/0","1|0")]=var_gr$REF_plus[var_gr$GT %in% c("1/0","1|0")]
  var_gr$genome2_minus[var_gr$GT %in% c("1/0","1|0")]=var_gr$REF_minus[var_gr$GT %in% c("1/0","1|0")]
  var_gr$genome2_tri[var_gr$GT %in% c("1/0","1|0")]=var_gr$REF_tri[var_gr$GT %in% c("1/0","1|0")]
  return(var_gr)
}

# Function to set genome and chromosome lengths to a GR object
setGenomeLengths <- function(GR,chrsOfInterest=paste("chr",1:22,sep=""),genome_in="hg19"){
  # Get genome info
  GR=chr_check(GR)
  genome.seqinfo <- seqinfo(getBSgenome(genome_in))
  genome.seqinfo <- genome.seqinfo[chrsOfInterest]
  GR <- GR[seqnames(GR) %in% chrsOfInterest]
  genome(GR) <- genome(genome.seqinfo)
  seqlevels(GR) <- seqlevels(genome.seqinfo)
  seqlengths(GR) <- seqlengths(genome.seqinfo)
  
  return(GR)
}
#read in all sample tissue diff
import.subject<-function(inDir,calc='diff'){
  #for calc: diff -> dMML etc, allele -> NME etc
  bed_in=dir(inDir,pattern="bedGraph")
  sample_in=unique(sub('_phased.*','',bed_in))
  GRs=GRanges()
  for (sp in sample_in) {
    subjects=sub('_.*','',sp)
    tissues=sub(paste0(subjects,"_"),'',sp)
    # Print sample being loaded
    cat("Loading:",subjects,'-',tissues,'\n')
    if (calc=='diff'){
      GR.in=read.diffGR(subjects,tissues,inDir)
    }else if(calc=='allele'){
      GR.in=read.alleleGR(subjects,tissues,inDir)
    }else {cat('Wrong calc \n')}
    GRs=c(GRs,GR.in)
  }
  return(GRs)
  
}
#Function to read in single GR object:
read.diffGR<-function(subjects,tissues,inDir,chrsOfInterest=paste("chr",1:22,sep="")){
  #Initialization
  GRs=GRanges()
  #Make sure the inputs are unique
  # dmml
  filename=paste(inDir,subjects,"_",tissues,"_phased_tmml_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=TRUE,'dMML')
  
  GRs <- c(GRs,GR.in)
  # dnme
  filename=paste(inDir,subjects,"_",tissues,"_phased_tnme_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=TRUE,'dNME')
  GRs <- c(GRs,GR.in)
  # uc
  filename=paste(inDir,subjects,"_",tissues,"_phased_tpdm_pvals.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=TRUE,'UC')
  GRs <- c(GRs,GR.in)
  #Check if pvalue available
  GRs <- GRs[!is.na(GRs$pvalue)]
  # Add sample field
  GRs$Sample <- paste(tissues,"-",GRs$Subject)
  # Add genome info 
  GRs <- setGenomeLengths(GRs,chrsOfInterest=chrsOfInterest)
  return(GRs)
}
#file_ends can be 
#c('dmml_pvals,dnme_pvals,uc_pvals')
#c('mml1,mml2,nme1,nme2')
read.alleleGR<-function(subjects,tissues,inDir,chrsOfInterest=paste("chr",1:22,sep="")){
  #MML1
  GRs=GRanges()
  filename=paste(inDir,subjects,"_",tissues,"_phased_mml1.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'MML')
  GR.in$Genome="1"
  GRs <- c(GRs, GR.in)
  #MML2
  filename=paste(inDir,subjects,"_",tissues,"_phased_mml2.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'MML')
  GR.in$Genome="2"
  GRs <- c(GRs, GR.in)
  #NME1
  filename=paste(inDir,subjects,"_",tissues,"_phased_nme1.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'NME')
  GR.in$Genome="1"
  GRs <- c(GRs, GR.in)
  #NME2
  filename=paste(inDir,subjects,"_",tissues,"_phased_nme2.bedGraph",sep="")
  GR.in=import.ASMbed(subjects,tissues,filename,pvalue=FALSE,'NME')
  GR.in$Genome="2"
  GRs <- c(GRs, GR.in)
  # Add sample field
  GRs$Sample <- paste(tissues,"-",subjects)
  # Add genome info 
  GRs <- setGenomeLengths(GRs,chrsOfInterest=chrsOfInterest)
  GRs$K=GRs$'NA.1'
  GRs$'NA.1'=NULL
  return(GRs)
}
#Read in each bed file, for new method, no need to resize
import.ASMbed<-function(subjects,tissues,filename,pvalue=TRUE,Statistic,chrsOfInterest=paste("chr",1:22,sep="")){
  GR <- import.bedGraph(filename)
  #fit  bedGraph reads, import.bedGraph will remove 1 from start
  #Check if all files are 0 based or 1 based? Check on genome browser, UCSC: check SNP location (0 based)
  start(GR)=start(GR)-1
  GR_out=chr_check(GR)
  GR_out$ASM=NULL
  GR_out$Data=NULL
  GR_out$Subject <- subjects
  GR_out$tissue<-tissues
  GR_out$Statistic <- Statistic
  GR_out$Value<- GR$score
  #for differential analysis, make sure the GR_out are unique
  if(pvalue){
    GR_out$pvalue <- as.numeric(GR$NA.)
    GR_out <- GR_out[!duplicated(GR[,c()])]
    }
  else{GR_out$N <- GR$NA.} #not diff for new samples
  return(GR_out)
}
#Get features
getGeneralFeats_CpG <- function(CpGdir,enhancerDir='',chrsOfInterest=paste("chr",1:22,sep="")){
  
  # Features included
  featureNickNames <- c("genome-wide","CpG island","CpG shore","CpG shelf","CpG open sea",
                        "gene body","exon","intron","intergenic")
  
  # Define list of feature GRs
  outGR <- GRangesList()
  GRtemp <- unlist(tileGenome(seqinfo(Hsapiens),ntile=1))
  
  outGR[["genome-wide"]] <- setGenomeLengths(GRtemp)
  #Redefining CpG islands using hidden Markov models 
  CpG_all <- readRDS(paste(CpGdir,"CpG_hg19.rds",sep=""))
  CpG_all<-setGenomeLengths(CpG_all)
  #Could also use UCSC genome browser CpG file
  cpg_islands <- readRDS(paste(CpGdir,"cpg_islands_hg19.rds",sep=""))
  cpg_islands<-subsetByOverlaps(CpG_all,cpg_islands)
  outGR[["CpG island"]] <- setGenomeLengths(cpg_islands)
  
  # extract the shore defined by 2000 bp upstream and downstream of cpg islands
  shore1 <- flank(cpg_islands, 2000)
  shore2 <- flank(cpg_islands,2000,FALSE)
  shore1_2 <- reduce(c(shore1,shore2))
  
  # extract the features (ranges) that are present in shores only and not in
  # cpg_islands (ie., shores not overlapping islands)
  cpgi_shores <- setdiff(shore1_2, cpg_islands)
  olap=findOverlaps(CpG_all,cpgi_shores)
  cpgi_shores<-subsetByOverlaps(CpG_all,cpgi_shores)
  outGR[["CpG shore"]] <- setGenomeLengths(cpgi_shores)
  
  # extract the shore defined by 4000 bp upstream and downstream of cpg islands
  shelves1 <- flank(cpg_islands, 4000)
  shelves2 <- flank(cpg_islands,4000,FALSE)
  shelves1_2 <- reduce(c(shelves1,shelves2))
  
  # create a set of ranges consisting CpG Islands, Shores
  island_shores <- c(cpg_islands,cpgi_shores)
  
  # extract the features (ranges) that are present in shelves only
  # and not in cpg_islands  or shores(ie., shelves not overlapping islands or shores)
  cpgi_shelves <- setdiff(shelves1_2, island_shores)
  cpgi_shelves<-subsetByOverlaps(CpG_all,cpgi_shelves)
  outGR[["CpG shelf"]] <- setGenomeLengths(cpgi_shelves)
  
  # Open sea
  open_sea <- setdiff(outGR[["genome-wide"]],c(outGR[["CpG island"]],outGR[["CpG shore"]],outGR[["CpG shelf"]]))
  open_sea<-subsetByOverlaps(CpG_all,open_sea)
  outGR[["CpG open sea"]] <- setGenomeLengths(open_sea)
  
  # Enhancers 
  #enhancers <- import.bed(paste(enhancerDir,"enhancers.bed",sep=""))[,c()]
  
  #outGR[["enhancer"]] <- setGenomeLengths(enhancers)
  
  # Other generic features
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes <- GenomicFeatures::genes(txdb)
  outGR[["gene body"]] <- setGenomeLengths(genes)
  exons <- GenomicFeatures::exons(txdb)
  outGR[["exon"]] <- setGenomeLengths(exons[,c()])
  introns <- GenomicFeatures::intronicParts(txdb)
  outGR[["intron"]] <- setGenomeLengths(introns[,c()])
  intergenic <- setdiff(outGR[["genome-wide"]],outGR[["gene body"]],ignore.strand=TRUE)
  outGR[["intergenic"]] <- setGenomeLengths(intergenic)
  #Use annotation hub for TSS, promoter have something to do with strand
  proms <- promoters(genes,upstream=2000,downstream=1000)#ask Michael about this
  outGR[["promoter"]] <- setGenomeLengths(proms)
  TSS<-promoters(genes,upstream=0,downstream=0)
  outGR[["TSS"]] <- setGenomeLengths(TSS)
  # Gene name mapping
  geneBodyNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["gene body"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["gene body"]]$gene_name <- geneBodyNameMap$SYMBOL
  promNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["promoter"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["promoter"]]$gene_name <- promNameMap$SYMBOL
  promNameMap <- AnnotationDbi::select(Homo.sapiens,key=as.character(outGR[["TSS"]]$gene_id),keytype="ENTREZID",columns=c("SYMBOL"))
  outGR[["TSS"]]$gene_name <- promNameMap$SYMBOL
  # Return
  return(outGR)
  
}

#Put hetCpG count into each sample in gr_allele

##redo instead of adding allele information to allele CpG, add it to GR_merge
#This is mainly modifed to fit the output of CPEL ASM

stat_merge<-function(gr_in,allele_in,vcf_in,CpG){
  #Check merge behavior, 
  dMML=gr_in[gr_in$Statistic=="dMML"]
  dNME=gr_in[gr_in$Statistic=="dNME"]
  UC=gr_in[gr_in$Statistic=="UC"]
  gr=unique(granges(gr_in))
  olap_dMML=findOverlaps(gr,dMML,type="equal")
  gr$dMML[queryHits(olap_dMML)]=dMML$Value[subjectHits(olap_dMML)]
  gr$dMML_pval[queryHits(olap_dMML)]=dMML$pvalue[subjectHits(olap_dMML)]
  gr$Sample=unique(gr_in$Sample)
  
  olap_dNME=findOverlaps(gr,dNME,type="equal")
  gr$dNME[queryHits(olap_dNME)]=dNME$Value[subjectHits(olap_dNME)]
  gr$dNME_pval[queryHits(olap_dNME)]=dNME$pvalue[subjectHits(olap_dNME)]
  
  olap_UC=findOverlaps(gr,UC,type="equal")
  gr$UC[queryHits(olap_UC)]=UC$Value[subjectHits(olap_UC)]
  gr$UC_pval[queryHits(olap_UC)]=UC$pvalue[subjectHits(olap_UC)]
  
  gr$NME1=gr$MML1=gr$NME2=gr$MML2=gr$N=NA
  #genome 1 NME
  olap_NME1=findOverlaps(gr,allele_in[allele_in$Genome==1 & allele_in$Statistic=='NME'],type="equal")
  gr$NME1[queryHits(olap_NME1)]=allele_in[allele_in$Genome==1 & allele_in$Statistic=='NME']$Value[subjectHits(olap_NME1)]
  gr$N[queryHits(olap_NME1)]=allele_in[allele_in$Genome==1 & allele_in$Statistic=='NME']$N[subjectHits(olap_NME1)]
  #genome 2 NME
  olap_NME2=findOverlaps(gr,allele_in[allele_in$Genome==2 & allele_in$Statistic=='NME'],type="equal")
  gr$NME2[queryHits(olap_NME2)]=allele_in[allele_in$Genome==2 & allele_in$Statistic=='NME']$Value[subjectHits(olap_NME2)]
  #genome 1 MML
  olap_MML1=findOverlaps(gr,allele_in[allele_in$Genome==1 & allele_in$Statistic=='MML'],type="equal")
  gr$MML1[queryHits(olap_MML1)]=allele_in[allele_in$Genome==1 & allele_in$Statistic=='MML']$Value[subjectHits(olap_MML1)]
  #genome 1 MML
  olap_MML2=findOverlaps(gr,allele_in[allele_in$Genome==2 & allele_in$Statistic=='MML'],type="equal")
  gr$MML2[queryHits(olap_MML2)]=allele_in[allele_in$Genome==2 & allele_in$Statistic=='MML']$Value[subjectHits(olap_MML2)]
  gr$Subject=unique(gr_in$Subject)
  gr$tissue=unique(gr_in$tissue)
  #Add g1cg etc information to the gr
  gr=add_hetCPG(gr,vcf_in,CpG)
  return(gr)
}

add_hetCPG<-function(gr,vcf_in,CpG){
  olap=findOverlaps(vcf_in,gr,type='within',select='all')
  #Count number of CG, here g1CG=genome1 and g2CG=genome2, however, we need to split based on the GT, also ref
  vcf_in$g1CG=as.numeric(grepl("CG",vcf_in$genome1_tri))
  vcf_in$g2CG=as.numeric(grepl("CG",vcf_in$genome2_tri))
  vcf_in$refCG=as.numeric(grepl("CG",vcf_in$REF_tri))
  vcf_in$altCG=as.numeric(grepl("CG",vcf_in$ALT_tri))
  #Count number of Het CpG here *Check 
  df_sub=data.table(subjHits=subjectHits(olap),
                    g1CG=vcf_in$g1CG[queryHits(olap)],g2CG=vcf_in$g2CG[queryHits(olap)],
                    refCG=vcf_in$refCG[queryHits(olap)],altCG=vcf_in$altCG[queryHits(olap)])
  
  agg_sub=df_sub[,.(g1CG=sum(g1CG),g2CG=sum(g2CG),refCG=sum(refCG),altCG=sum(altCG)),by=subjHits]
  gr$g1CG[agg_sub$subjHits]=agg_sub$g1CG#agg_sub$subjHits is the unique subject hits
  gr$g2CG[agg_sub$subjHits]=agg_sub$g2CG
  gr$refCG[agg_sub$subjHits]=agg_sub$refCG#agg_sub$subjHits is the unique subject hits
  gr$altCG[agg_sub$subjHits]=agg_sub$altCG
  gr$N_hg19=countOverlaps(gr,CpG)
  #count number of CG lost from reference
  gr$N_nonhet=gr$N_hg19-countOverlaps(gr,vcf_in[(vcf_in$refCG-vcf_in$altCG)>0])
  
  return(gr)
}


#add genes to GR
add_gene_GR<-function(GR_merge_in,feature_in,feature_names){
  mcols(GR_merge_in)[[feature_names]]=NA
  olap_promoter=findOverlaps(GR_merge_in,feature_in)
  df_idx=data.table(qt=queryHits(olap_promoter),
                    genes=as.character(feature_in$gene_name[subjectHits(olap_promoter)]),
                    stringsAsFactors = F)
  df_idx=df_idx[,list(genes=list(genes)),by=list(df_idx$qt)]
  colnames(df_idx)=c('qt','genes')
  df_idx$genes=lapply(df_idx$genes,as.character)
  mcols(GR_merge_in)[[feature_names]][df_idx$qt]=df_idx$genes
  return(GR_merge_in)
  
}

#Count Number of Het CpG at extened regions, each extend 500 bp
hetCGallele_merged<-function(sub,gr_merge,CpG,vcf_in,gene_size=500){
  cat('Analyzing',sub,'\n')
  t1=proc.time()[[3]]
  #Import vcf file
  sub_vcf=vcf_in[[sub]]
  sub_allele=gr_merge[gr_merge$Subject==sub]
  #gr_allele got resized with start and end +1, use +2 to include equal start & end, no longer needed for new output?
  #sub_het=resize(sub_het, width(sub_het) + 4, fix="center")
  sub_allele$gff_size=width(sub_allele)
  gr_out=GR_resize_merged(sub_allele,CpG,sub_vcf,gene_size=gene_size,sub)#Check for calculating density
  cat('Finish analyzing',sub,proc.time()[[3]]-t1,'\n')
  return(gr_out)
}
#Generate density from GR allele, at least need 200 bp for density
#Count number of hetCG at each allele GR from result, add gff size, and N
GR_resize_merged<-function(GR.in,CpG_sites,hetCpG,sub,gene_size=500){
  ##From definitino of CpG island, use 200 bp regions
  GR.extend=resize(GR.in,width=width(GR.in)+gene_size,fix='center')
  GR.in$CG_hg19_extend=countOverlaps(GR.extend,CpG_sites)
  #any SNP contain CG in ref genome
  GR.in$CG_nonhet_extend=GR.in$CG_hg19_extend-countOverlaps(GR.extend,hetCpG[grepl("CG",hetCpG$REF_tri)])
  gr_seq=getSeq(Hsapiens,GR.extend,as.character=T)
  GR.in$CGcont_exp=do.call('c',lapply(gr_seq,countCGOR))
  #Count CpG in genome 1
  GR.in$CG_allele_extend_g1=GR.in$CG_nonhet_extend+
    countOverlaps(GR.extend,hetCpG[grepl("CG",hetCpG$genome1_tri)])
  GR.in$CG_allele_extend_g2=GR.in$CG_nonhet_extend+
    countOverlaps(GR.extend,hetCpG[grepl("CG",hetCpG$genome2_tri)])
  GR.in$gff_size_extend=width(GR.extend)
  return(GR.in)
}

#Count expted CG numbers and return data.frame
countCGOR<-function(x){ #x=input seq
  #calculate Odds ratio for expected CG vs actual CG
  #Expected CG number C * number G/total length
  # Gardiner-Garden M, Frommer M (1987). "CpG islands in vertebrate genomes". Journal of Molecular Biology.
  #Wiki, actual: ((number of C + Number of G)/2)^2/length of genomics Normalized CpG content, whole genome ~25%
  NC=countPattern('C',x)
  NG=countPattern('G',x)
  #NCG=countPattern('CG',x)
  CG_exp=NC*NG/nchar(x) #PC*PG*length
  CG_exp_norm=((NC+NG)/2)^2/nchar(x) #assuming PC=PG = (NC+NG)/2/length
  return(CG_exp_norm)
}


#Give each SNP an ASM information for each subject
variant_meta<- function(subj,variant_in,GR_in){ #variant_in for each subject, GR_in for each subject
  cat('Processing',subj,'\n')
  GR_in_subj=GR_in[GR_in$Subject==subj]
  variant_in_subj=variant_in[[subj]]
  sp=unique(GR_in_subj$Sample)
  gr_out_sp = GRanges()
  for (sps in sp){
    cat('Processing',sps,'\n')
    gr_out_sp=c(gr_out_sp,variant_meta_sp(variant_in_subj,GR_in_subj[GR_in_subj$Sample==sps]))
    
  }
  return(gr_out_sp)
  
  
}

#For each sample, assign NME to SNP within each region
variant_meta_sp <-function(variant_subj,GR_sp){
  olap=findOverlaps(variant_subj,GR_sp,maxgap =0,type='within')
  gr_out=variant_subj[queryHits(olap)]
  #Find GR_merge olap
  olap=subjectHits(olap)
  mcols(gr_out)=cbind(mcols(gr_out),mcols(GR_sp)[olap,])
  gr_out$refNME=NA
  gr_out$altNME=NA
  gr_out$refMML=NA
  gr_out$altMML=NA
  #NME
  gr_out$refNME[gr_out$GT %in% c("0/1","0|1")]= gr_out$NME1[gr_out$GT %in% c("0/1","0|1")]
  gr_out$altNME[gr_out$GT %in% c("0/1","0|1")]= gr_out$NME2[gr_out$GT %in% c("0/1","0|1")]
  gr_out$refNME[gr_out$GT %in% c("1/0","1|0")]= gr_out$NME2[gr_out$GT %in% c("1/0","1|0")]
  gr_out$altNME[gr_out$GT %in% c("1/0","1|0")]= gr_out$NME1[gr_out$GT %in% c("1/0","1|0")]
  
  #MML
  gr_out$refMML[gr_out$GT %in% c("0/1","0|1")]= gr_out$MML1[gr_out$GT %in% c("0/1","0|1")]
  gr_out$altMML[gr_out$GT %in% c("0/1","0|1")]= gr_out$MML2[gr_out$GT %in% c("0/1","0|1")]
  gr_out$refMML[gr_out$GT %in% c("1/0","1|0")]= gr_out$MML2[gr_out$GT %in% c("1/0","1|0")]
  gr_out$altMML[gr_out$GT %in% c("1/0","1|0")]= gr_out$MML1[gr_out$GT %in% c("1/0","1|0")]
  
  gr_out$Statistic=GR_sp$Statistic[olap]
  
  # gr_out$HetCpG=FALSE
  # #This may be different from previous result?
  # gr_out$HetCpG=((gr_out$REF_plus =='CG' | gr_out$REF_minus=='CG') & !(gr_out$ALT_plus =='CG' | gr_out$ALT_minus=='CG')) |
  #   (!(gr_out$REF_plus =='CG' | gr_out$REF_minus=='CG') & (gr_out$ALT_plus =='CG' | gr_out$ALT_minus=='CG'))
  return(gr_out)
}

#Read in allele-agnositc model
read.agnostic<-function(file_in,GR_merge_in=NULL,allele_include=T,olap_type="any",all_regions=NA,sample_in=NA,hyper_var_file=NA){
  stat=toupper(sub('.*_','',sub('.bedGraph','',file_in)))
  informME_in=read.bedGraph.informME(file_in)
  if(length(GR_merge_in)>0){
    GR_merge_in=GR_merge_in[GR_merge_in$Sample==sample_in]
  #Find overlapped region
    olap=findOverlaps(informME_in,GR_merge_in,type=olap_type)
    if(length(olap)>0) {informME_in=informME_in[-queryHits(olap)]}
  }
  #add GR_merge data
  if(allele_include){
  
    cat('Percent overlap with dNME region:',length(unique(queryHits(olap)))/length(informME_in)*100,'%\n')
    informME_in$score_original=informME_in$score
 
    #replace value instead of remove regions
    olap=findOverlaps(all_regions,GR_merge_in)
    asm_replacement=data.table(idx=queryHits(olap), 
                               score= rowMeans(as.matrix(elementMetadata(GR_merge_in)[paste(stat,c('1','2'),sep='')]))[subjectHits(olap)],
                               dMML=GR_merge_in$dMML[subjectHits(olap)],
                               dMML_pval=GR_merge_in$dMML_pval[subjectHits(olap)],
                               N=GR_merge_in$N[subjectHits(olap)])
    asm_replacement=asm_replacement[,list(score=mean(score),dMML=mean(dMML),dMML_pval=mean(dMML_pval)),by=list(idx)]
    all_regions=all_regions[asm_replacement$idx]
    all_regions$score=asm_replacement$score
    all_regions$K=NA
    all_regions$dMML=asm_replacement$dMML
    all_regions$dMML_pval=asm_replacement$dMML_pval
    informME_in=c(informME_in,all_regions)
  }
 informME_in=informME_in[!is.infinite(informME_in$score)]
  informME_in$Sample=sample_in
  informME_in$hyper_var_fn=hyper_var_file
  return(informME_in)

}
#read in hypervaribility from file
read_hypervar<-function(hyper_var_in){
  hyper_var=data.table()
  for (fn in hyper_var_in){
    fn_in=readRDS(fn)
    fn_in$gene_name=rownames(fn_in)
    fn_in=as.data.table(fn_in)
    rownames(fn_in)=NULL
    hyper_var=rbind(hyper_var,fn_in)
  }
  
  #Find the dataset that have most regions overlap with NME regions, check correlation?
  hyper_var=hyper_var[,list(mean=mean(mean,na.rm=T),var=mean(var,na.rm=T),
                            hypervar_var=mean(hypervar_var,na.rm=T),hypervar_logvar=mean(hypervar_logvar,na.rm=T)),
                      by=list(hyper_var$gene_name)]
  colnames(hyper_var)[1]='gene_name'
  return(hyper_var)
}
#Calculating NME vs hypervaribility
dist_plot_calc<-function(informME_in,genes_hypervar,genomic_features,enhancer=FALSE){
  if(enhancer){informME_in_dist=informME_in}else{informME_in_dist=dist_calc(informME_in,genomic_features$TSS)}
  mcols(informME_in_dist)=cbind(mcols(informME_in_dist),genes_hypervar[,-1][match(informME_in_dist$gene,genes_hypervar$gene_name)])
  
  return(informME_in_dist)
}
#calculate distance to the genomic features
dist_calc<-function(gr_in,gr_feature){
  
  dist_nearest=nearest(resize(gr_in,1,fix="center"),resize(gr_feature,1,fix="start"))
  
  gr_non_na=which(!is.na(dist_nearest))
  gr_feature=gr_feature[dist_nearest[gr_non_na]]
  gr_in=gr_in[gr_non_na]
  gr_in$dist=NA
  gr_in$gene=NA
  
  sgn <- as.integer(ifelse(strand(gr_feature)=="+",1,-1))
  gr_in$dist=sgn*(start(resize(gr_in,1,fix="center"))-start(resize(gr_feature,1,fix="start")))
  gr_in$gene=gr_feature$gene_name
  
  return(gr_in)
}
read.bedGraph.informME<-function(file_in){
  
  informME_in=import.bedGraph(file_in)
  if(length(informME_in)>0){
    colnames(elementMetadata(informME_in))=c('score','N','K')
    if(all(seqlevels(informME_in)==gsub('chr','',seqlevels(informME_in)))){seqlevels(informME_in)=paste('chr',seqlevels(informME_in),sep='')}
    #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(informME_in)=start(informME_in)-1
    return(informME_in)
  }
}
read_chromHMM_bed_PRC<-function(bed_dir,rep,blacklist_region,cpgr){
  bed_out=GRanges()
  for(fn in dir(bed_dir,pattern='.bed.gz')){
    #get sample name etc
    fn_split=strsplit(fn,'_')[[1]]
    stage=gsub('e','day',fn_split[1])
    stage=gsub('P','day',stage)
    stage=gsub('\\.','\\_',stage)
    tissue=gsub('facial-prominence','EFP',fn_split[2])
    tissue=gsub('neural-tube','NT',tissue)
    bed_in=read.table(paste(bed_dir,fn,sep=''))
    colnames(bed_in)=c('chr','start','end','chrom_num','chrom_state')
    bed_in=makeGRangesFromDataFrame(bed_in,keep.extra.columns = T)
    bed_in=reduce(bed_in[bed_in$chrom_state%in%c('Hc-P','Pr-B')])
    bed_in$stage=stage
    bed_in$tissue=tissue
    bed_in$rep=rep
    bed_in$Sample=paste(stage,tissue,sep='-')
    bed_out=c(bed_out,bed_in)
    
  }
  bed_out=subdivideGRanges(bed_out,250)
  olap=findOverlaps(bed_out,cpgr)
  bed_out=unique(bed_out[queryHits(olap)])
  
  olap=findOverlaps(bed_out,blacklist_region)
  bed_out=bed_out[-queryHits(olap)]
  return(bed_out)
}
#read in mouse enhancer
read_chromHMM_bed<-function(bed_dir,rep){
  bed_out=GRanges()
  for(fn in dir(bed_dir,pattern='.bed.gz')){
    #get sample name etc
    fn_split=strsplit(fn,'_')[[1]]
    stage=gsub('e','E',fn_split[1])
    tissue=gsub('facial-prominence','EFP',fn_split[2])
    tissue=gsub('neural-tube','NT',tissue)
    bed_in=read.table(paste(bed_dir,fn,sep=''))
    colnames(bed_in)=c('chr','start','end','chrom_num','chrom_state')
    bed_in=makeGRangesFromDataFrame(bed_in,keep.extra.columns = T)
    bed_in$stage=stage
    bed_in$tissue=tissue
    bed_in$rep=rep
    bed_in$Sample=paste(stage,tissue,sep='-')
    bed_out=c(bed_out,bed_in)
  }
  return(bed_out)
}
#Enrichment analysis for genomic features in dMML and dNME
testEnrichmentFeature_stat<-function(dataGR,featureGR,maxgap=0,output_ft=1){
  # Find ranges overlapping with feature
  olaps <- findOverlaps(dataGR,featureGR,type="any",select="all",maxgap = maxgap)
  
  indFeature <- queryHits(olaps)
  featurestatistic <- dataGR[indFeature]
  complementarystatistic <- dataGR[-indFeature]
  
  # Enrichment of in feature
  #featurestatistic <- featureData[featureData$Statistic==statistic]
  #complementarystatistic <- complementaryData[complementaryData$Statistic==statistic]
  contTablestatistic <- data.frame(ASM=c(NA,NA),nonASM=c(NA,NA))
  rownames(contTablestatistic) <- c("Feature","Complementary")
  contTablestatistic[1,]$ASM <- sum(featurestatistic$ASM=="Yes")
  contTablestatistic[1,]$nonASM <- sum(featurestatistic$ASM=="No")
  contTablestatistic[2,]$ASM <- sum(complementarystatistic$ASM=="Yes")
  contTablestatistic[2,]$nonASM <- sum(complementarystatistic$ASM=="No")
  #print(contTablestatistic)
  #Return overlap >=output_ft
  if(contTablestatistic[1,1]>=output_ft){
  return(list(contTablestatistic,fisher.test(contTablestatistic)))
  }
}
#MAE enrich
MAE_enrich<-function(GR_merge,pval_cutoff,genes='genes_promoter',stat='dMML_pval',MAE=MAE){
  #GR_merge=GR_merge[!is.na(GR_merge$genes_promoter)]
  GR_merge=elementMetadata(GR_merge)
  stat_gene=GR_merge[[genes]][GR_merge[[stat]]<=pval_cutoff]
  non_stat_gene=GR_merge[[genes]][GR_merge[[stat]]>pval_cutoff]
  stat_MAE=sum(unlist(lapply(stat_gene,function(x) any(x %in% MAE))))
  nonstat_MAE=sum(unlist(lapply(non_stat_gene,function(x) any(x %in% MAE)))) 
  nonstat_nonMAE=sum(!unlist(lapply(non_stat_gene,function(x) any(x %in% MAE))))
  stat_nonMAE=sum(!unlist(lapply(stat_gene,function(x) any(x %in% MAE))))
  print(matrix(c(stat_MAE,stat_nonMAE,nonstat_MAE,nonstat_nonMAE),nrow=2))
  ft=fisher.test(matrix(c(stat_MAE,stat_nonMAE,nonstat_MAE,nonstat_nonMAE),nrow=2))
  return(data.frame(OR=ft$estimate,pvalue=ft$p.value,lowerCI= ft$conf.int[1],upperCI=ft$conf.int[2]))
}


#ChromHMM* check
ENCODE_to_sample<-function(sample_in){
  #Assign code to each sample
  ENCODE_number_to_sample=data.frame(sample=sample_in,stringsAsFactors = F)
  ENCODE_number_to_sample$subject=unlist(lapply(ENCODE_number_to_sample$sample,function(x) strsplit(x,' - ')[[1]][2]))
  ENCODE_number_to_sample$tissue=unlist(lapply(ENCODE_number_to_sample$sample,function(x) gsub(c('_single','_paired'),'',strsplit(x,' - ')[[1]][1])))
  #Assign known code to given sample
  ENCODE_number_to_sample$ENCODE=NA
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$subject=='H1']='E003'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='stem_27_undifferentiated_paired - HUES64']='E016'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='42_embryonic_stem_cell_single - H9']='E008'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='foreskin_melanocyte_paired - skin03']='E061'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='foreskin_keratinocyte_paired - skin03']='E058'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Small_Intestine']='E109'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Gastric']='E094'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Left_Ventricle']='E095'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Lung']='E096'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Psoas_Muscle']='E100'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Right_Ventricle']='E105'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Sigmoid_Colon']='E106'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Spleen']='E113'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Thymus']='E112'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue %in% c('Adipose','Adipose_Tissue')]='E063'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Aorta']='E065'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Esophagus']='E079'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Adrenal_Gland']='E080'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Ovary']='E097'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Pancreas']='E087'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Liver']='E066'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$tissue=='Right_Atrium']='E014'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='ectoderm_paired - HUES64']='E012'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='endoerm_27_paired - HUES64']='E011'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='mesoderm_23_paired - HUES64']='E013'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='brain_germinal_matrix_tissue_paired - HuFGM02']='E070'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='Brain_substantia_nigra_paired - 112']='E074'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='Brain_Hippocampus_middle_paired - 149']='E071'
  ENCODE_number_to_sample$ENCODE[ENCODE_number_to_sample$sample=='Brain_Hippocampus_middle_paired - 150']='E071'
  
  return(ENCODE_number_to_sample)
}
##Use  CMH test for chromHMM analysis
chromHMM_OR<-function(GR_merge,chromHMM,sample_name,pval_cutoff=0.1,stat="dNME_pval"){
  
  GR_merge_sp=GR_merge[GR_merge$Sample%in%sample_name]
  GR_merge_sp$ASM="No"
  GR_merge_sp$ASM[elementMetadata(GR_merge_sp)[,stat]<=pval_cutoff]="Yes"
  out_df=data.frame()
  count_table=list()
  count_table_N=data.frame()
  for(states in unique(chromHMM$name)){
    
    OR=testEnrichmentFeature_stat(GR_merge_sp,chromHMM[chromHMM$name==states])
    if(!is.null(OR)){
    #Get contengency table
    count_table[[states]]=OR[[1]]
    OR=OR[[2]]
    #result directly from fisher.test for each state
    out_df=rbind(out_df,
                 data.frame(state=states,OR=OR$estimate,p_value=OR$p.value,
                            lower_CI=OR$conf.int[1],upper_CI=OR$conf.int[2]))
    }
    # }
  }
  return(list(out_df,count_table))
}
#Combine contengency table
chromHMM_combine<-function(chromHMM_in){
  cont_table_all=list()
  #ChromHMM, list object, each is a sample, within sample ,each is a state
  for (states in unique(unlist(lapply(chromHMM_dMML_all_ls,function(x) names(x[[2]]))))){
    #extract 2x2 table for each states for each sample, return a list of sample with its contengency table
    chromHMM_in_cont=lapply(chromHMM_in, function(x,states) x[[2]][[states]],states=states)
    #For each state, construct a CMH table, 4 columns
    cont_table_all[[states]]=do.call(rbind,lapply(seq_along(chromHMM_in_cont),melt_cont,cont_in=chromHMM_in_cont))
    print(states)
    
  }
  cont_table_all_CMH=lapply(cont_table_all,CMH_test)
  cont_table_all_CMH_df=data.frame(states=names(cont_table_all_CMH), 
                                   OR=unlist(lapply(cont_table_all_CMH, function(x) x$estimate[[1]])), 
                                   p_value=unlist(lapply(cont_table_all_CMH, function(x) x$p.value)),
                                   lower_CI=unlist(lapply(cont_table_all_CMH,function(x) x$conf.int[1])),
                                   upper_CI=unlist(lapply(cont_table_all_CMH,function(x) x$conf.int[2])))
  
  return(cont_table_all_CMH_df)
}

melt_cont<-function(i,cont_in){
  sp=names(cont_in)[i]
  # cont_out=cont_in[[i]]
  # cont_out$subject=factor(paste(sp,cont_out$N,sep='_'))
  # return(cont_out)
  count_out=as.numeric(unlist(cont_in[[i]]))
  if(length(count_out)>0){
  #print(count_out)
  #cutoff of smallest number in the table
  if (all(count_out>0)){
    return(data.frame(subject=factor(rep(sp,4)),
                      ASM=factor(c('ASM','ASM','Non-ASM','Non-ASM')),
                      feature=factor(c('Feature','Non_feature','Feature','Non_feature')),
                      count=count_out
    ))
  }
  }
}

CMH_test<-function(df_in,CMH_eqn=count~ASM+feature+subject){
  if (nrow(df_in)>0){
    CMH_table=xtabs(CMH_eqn,data=df_in)
    
    if (length(dim(CMH_table)[3]!=0)){
      
      if(dim(CMH_table)[3]>1){
        fs=mantelhaen.test(CMH_table)
      }
      else if(dim(CMH_table)[3]==1){
        #print(as.matrix(CMH_table[,,1]))
        #print(count_all)
        fs=fisher.test(as.matrix(CMH_table[,,1]))
        
      }
      return(fs)
    }
    
  }
}
#SNP OR
OR_calc<-function(tb_in,SNP,SNP_name,pval_cutoff=NA,stat_in="NME"){
 print(SNP)
  larger_SNP=sum(tb_in[[SNP_name]]==SNP&(tb_in[[paste0('d',stat_in,"_relative")]]>0))
  larger_nonSNP=sum(tb_in[[SNP_name]]!=SNP&(tb_in[[paste0('d',stat_in,"_relative")]]>0))
  lower_SNP=sum(tb_in[[SNP_name]]==SNP&!(tb_in[[paste0('d',stat_in,"_relative")]]>0))
  lower_nonSNP=sum(tb_in[[SNP_name]]!=SNP&!(tb_in[[paste0('d',stat_in,"_relative")]]>0))
  cont_table=matrix(c(larger_SNP,larger_nonSNP,lower_SNP,lower_nonSNP),nrow=2)
  print(cont_table)
  OR=fisher.test(cont_table,alternative ='two.sided')
  #print(cont_table)
  return(data.table(OR=OR$estimate,pvalue=OR$p.value,lowerCI=OR$conf.int[1],upperCI=OR$conf.int[2],SNP=SNP))
}

dist_plot_run<-function(informME_in_dist,theme_glob,ylab,stat_in,cutoff=pval_cutoff,dir){
  
  #informME_in_dist=informME_in_dist[-which(informME_in_dist$dMML_pval<=cutoff)]
  
  informME_in_dist$exp_stat=informME_in_dist[[stat_in]]
  informME_in_dist=informME_in_dist[!is.na(exp_stat)]
  plot_informME_dat=informME_in_dist[,list(score=score,stat_in=stat_in,dist=dist,exp_stat=exp_stat,gene=gene,
                                           quant=findInterval(exp_stat,quantile(unique(data.table(gene=gene,exp_stat=exp_stat))$exp_stat,prob=c(0,0.25,0.5,0.75),na.rm=T)),
                                           hypervarquant=findInterval(exp_stat,quantile(unique(data.table(gene=gene,exp_stat=exp_stat))$exp_stat,prob=seq(0.01,1,0.01),na.rm=T))/100),
                                     #scorequant001=findInterval(exp_stat,quantile(score,prob=c(0.01,1,0.01),na.rm=T))/100
                                     by=list(Sample)]
  rm(informME_in_dist)
  quant=c("0-25%","25%-50%","50%-75%","75%-100%")
  plot_informME_dat$quant=quant[plot_informME_dat$quant]
  plot_informME_dat$quant=as.factor(plot_informME_dat$quant)
  #dist_plot=list()t
  #print(unique(plot_informME_dat[,c(3,4)]))
  pdf(paste0(dir,'Figure3A_',ylab,'_',stat_in,'_dist.pdf'),width=3.5,height=3.5)
  for(sp in unique(plot_informME_dat$Sample)){
  print(ggplot(plot_informME_dat[abs(dist)<=3000&Sample==sp],aes(x=dist,y=score,color=quant))+theme_glob+
    geom_smooth(size=1,na.rm=TRUE,se=TRUE)+ggtitle(sp)+
    theme(legend.position="bottom")+ labs(color = "quantile")+
    xlab("Distance to TSS")+ylab(ylab)+guides(color=guide_legend(nrow=2,byrow=TRUE)))
  }
 dev.off()
  #find where min NME is 
  #plot_informME_dat_median=plot_informME_dat[abs(dist)<=3000,list(median_score=median(score)),by=list(dist,quant)]
  #min_dist=plot_informME_dat_median[which.min(median_score)]$dist
  #plot_informME_dat_ft=plot_informME_dat[dist<=min_dist+250&dist>=min_dist-250]
  plot_informME_dat_ft=plot_informME_dat[dist>=0&dist<=500]
  #print(cor.test(plot_informME_dat_cor$score,plot_informME_dat_cor$exp_stat))
  #Use correlation between hypervaribility quantile and NME
  # plot_informME_dat_cor=plot_informME_dat_ft[,list(median_score=median(score)),
  #                                         by=list(exp_stat,gene,hypervarquant,Sample)]

  print(plot_informME_dat_ft)
  cor_mean=plot_informME_dat_ft[,list(cor=cor(score,hypervarquant),
                                      pvalue=cor.test(score,hypervarquant)$p.value,
                                      lowerCI=cor.test(score,hypervarquant)$conf.int[1],
                                      upperCI=cor.test(score,hypervarquant)$conf.int[2],
                                      df=cor.test(score,hypervarquant)$parameter,
                                      stat=cor.test(score,hypervarquant)$statistic),by=list(Sample)]
  write.csv(cor_mean,paste0(dir,'Figure3A_',ylab,'_',stat_in,'_500bp_tb.csv'),row.names = F)
  print(mean(cor_mean$cor))
  # print(cor.test(plot_informME_dat_cor$median_score,plot_informME_dat_cor$hypervarquant))
  heatmap_informME_dat=plot_informME_dat_ft[,list(median_score=median(score)),by=list(hypervarquant,Sample)]
  heatmap_informME_dat=heatmap_informME_dat[,list(hypervarquant=hypervarquant,NME=median_score,cor=cor(median_score,hypervarquant)),by=list(Sample)]
  heatmap_informME_dat=heatmap_informME_dat[order(heatmap_informME_dat$cor,decreasing = F),]
  heatmap_informME_dat$Sample=factor(heatmap_informME_dat$Sample,levels = unique(heatmap_informME_dat$Sample))
  heatmap_plot=ggplot(heatmap_informME_dat,aes(hypervarquant,Sample,fill=NME))+geom_tile()+scale_fill_distiller(palette = "RdPu", direction = 1)+
    xlab('quantile')+ylab('Sample')+theme_glob+theme(legend.position = 'bottom')
  ggsave(paste0(dir,'Figure3A_',ylab,'_',stat_in,'_heatmap.pdf'),heatmap_plot,device='pdf',width=7,height=7,units="in")
  # cor_plot=ggplot(plot_informME_dat,aes(x=exp_stat,y=dat))+
  #   geom_smooth(size=1,na.rm=TRUE,se=TRUE)+theme_glob+
  #   theme(legend.position="bottom")+ggtitle('Correlation to hypervaribility within 500 bp upstream from TSS')+geom_point(alpha=0.1)
  # 
  # print(cor_plot)
  return(plot_informME_dat)
}
direction_calc_enriched_subj<-function(motif_loc,variant_all,gene_in,pval_cutoff=0.1,stat="NME"){
  
  #print(motif_sig_df)
  # motif_sig_df=motif_sig_df[which(motif_sig_df$qval<=0.2&motif_sig_df$OR>1),]
  motif_direction_out=mclapply(gene_in,direction_enriched_sample,
                               variant_gene=variant_all,motif_gene_subj=motif_loc,pval_cutoff=pval_cutoff,stat=stat,mc.cores =10)
  
  
  
  return(do.call(rbind,motif_direction_out))
}
merge_SNP_motif<-function(variant_gene,motif_gene_subj,stat){
  olap=findOverlaps(variant_gene,motif_gene_subj)
  variant_gene=variant_gene[queryHits(olap)]
  variant_gene$alleleDiff=motif_gene_subj$alleleDiff[subjectHits(olap)]
  variant_gene=variant_gene[!is.na(variant_gene$alleleDiff)]
  variant_gene$stat_diff=mcols(variant_gene)[[paste0('alt',stat)]]-mcols(variant_gene)[[paste0('ref',stat)]]

  return(variant_gene)
  
}
plot_merge_SNP_motif<-function(variant_gene,motif_gene_subj,stat,motif,pval_cutoff,theme_glob){
  motif_out=merge_SNP_motif(variant_gene[mcols(variant_gene)[[paste0('d',stat,'_pval')]]<=pval_cutoff],motif_gene [motif_gene$geneSymbol==motif],stat=stat)
  motif_out=convert_GR(motif_out,direction='DT')
  motif_out$`High binding affinity`=NaN
  motif_out$`Low binding affinity`=NaN
  #AlleleDiff is using alt-ref
  motif_out[alleleDiff<0]$`High binding affinity`=motif_out[alleleDiff<0,paste0("ref",stat),with=F]
  motif_out[alleleDiff>0]$`High binding affinity`=motif_out[alleleDiff>0,paste0("alt",stat),with=F]
  motif_out[alleleDiff<0]$`Low binding affinity`=motif_out[alleleDiff<0,paste0("alt",stat),with=F]
  motif_out[alleleDiff>0]$`Low binding affinity`=motif_out[alleleDiff>0,paste0("ref",stat),with=F]
  motif_out=melt.data.table(motif_out[,list(`High binding affinity`,`Low binding affinity`,region)],id.vars='region',value.name = "stat")
  print(ggplot(motif_out,aes(x=variable,y=stat))+
          geom_violin()+xlab("")+
          stat_summary(fun=median, geom="point", size=2, color="red")+
          ylab(stat)+ggtitle(motif)+theme_glob)
  return(motif_out)
}
direction_enriched_sample<-function(tf,variant_gene,motif_gene_subj,pval_cutoff,nperm=0,stat="NME"){
  variant_gene=merge_SNP_motif(variant_gene[mcols(variant_gene)[[paste0('d',stat,'_pval')]]<=pval_cutoff],
                                motif_gene_subj[motif_gene_subj$geneSymbol==tf],stat=stat)
  #alleleDiff is calculated use ref - alt, prefer low ent ones

  same_dir=sum(sign(variant_gene$alleleDiff)== sign(variant_gene$stat_diff),na.rm = TRUE)
  opposite_dir=sum(sign(variant_gene$alleleDiff)!= sign(variant_gene$stat_diff),na.rm = TRUE)
  # same_dir=sum(sign(variant_gene$alleleDiff)== sign(variant_gene$MML_diff),na.rm = TRUE)
  # opposite_dir=sum(sign(variant_gene$alleleDiff)!= sign(variant_gene$MML_diff),na.rm = TRUE)
  total_data=same_dir+opposite_dir
  
  variant_gene_df=data.frame(alleleDiff=sign(variant_gene$alleleDiff),stat_diff=sign(variant_gene$stat_diff))
  len_x=nrow(variant_gene_df)
  if(nperm>0){
    same_dir_perm=replicate(nperm,
                            sum(sample(variant_gene_df$alleleDiff,len_x,replace = F)==sample(variant_gene_df$stat_diff,len_x,replace = F)))
    same_dir_perm_prob=same_dir_perm/total_data
  }else{same_dir_perm_prob=-1}
  if(same_dir >0 |opposite_dir>0){
    
    binom=binom.test(same_dir,(same_dir+opposite_dir),0.5)
    #print(binom)
    # binom=summary(lm(abs(variant_gene$alleleDiff)~abs(variant_gene$dMML)))
    #binom$p.value=pf(binom$fstatistic[1],df1 = binom$fstatistic[2],df2 = binom$fstatistic[3],lower.tail = F)
    # return(data.frame(TF=unique(motif_gene_subj$geneSymbol),total_data=same_dir+opposite_dir,same_dir=same_dir,opposite_dir=opposite_dir,
    #                   binom.pval=binom$p.value,prob=binom$estimate[[1]],NSNP=length(variant_gene),stringsAsFactors = F))
    prob_binom=binom$estimate[[1]]
    if(nperm>0){binom.pval=sum(abs(same_dir_perm_prob-0.5)>=abs(prob_binom-0.5))/nperm}else{binom.pval=NA}
    # if(prob_binom>0.5){
    #   binom.pval=sum(same_dir_perm_prob>=(same_dir/total_data)|same_dir_perm_prob<=(1-(same_dir/total_data)))/nperm
    #   
    # }else if(prob_binom<0.5){binom.pval=sum(same_dir_perm_prob<=(same_dir/total_data)|same_dir_perm_prob>=(1-(same_dir/total_data)))/nperm}
    #cat(tf,':',binom.pval,'\n')
    return(data.table(TF=tf,total_data=total_data,same_dir=same_dir,opposite_dir=opposite_dir,
                      binom.pval_perm=binom.pval,binom.pval=binom$p.value,prob=prob_binom,NSNP=length(variant_gene),
                      lowerCI=binom$conf.int[1],upperCI=binom$conf.int[2],stringsAsFactors = F))
  }
}
OMIM_annotation<-function(motif_in,OMIM){
  motif_in=motif_in[order(motif_in$Proportion,decreasing=T),]
  motif_in$OMIM=NA
  for(tf in unique(motif_in$TF)){
    tf_in=gsub('\\(var.2\\)','',tf)
    tf_in=gsub('\\(var.3\\)','',tf_in)
    
    tf_in=unlist(strsplit(tf_in,'::'))
    #print(tf_in)
    OMIM_disease=OMIM$Phenotypes[which(unlist(lapply(OMIM$`Gene Symbols`,function(x) any(x%in% tf_in))))]
    if(length(OMIM_disease)>0){
      motif_in$OMIM[motif_in$TF==tf]=paste(OMIM_disease,collapse=';')}
    
  }
  
  return(motif_in)
}
NME_dNME_ken<-function(motif_in,GR_in,stat_in){
  tt1=proc.time()[[3]]
  olap=findOverlaps(motif_in,GR_in,select='all')
  if(length(GR_in$NME1)>0&length(GR_in$MML1)){
    GR_in$NME=(GR_in$NME1+GR_in$NME2)/2
    GR_in$MML=(GR_in$MML1+GR_in$MML2)/2
    }
  subj_olap=subjectHits(olap)
  stat_in_df=data.table(qt=queryHits(olap),stat_in=elementMetadata(GR_in)[[stat_in]][subj_olap],
                        sample=GR_in$Sample[subj_olap],stringsAsFactors = T)
  sample_all=unique(GR_in$Sample)
  gc()
  stat_in_df_stat=dcast.data.table(stat_in_df,qt~sample,value.var = "stat_in",
                                   fun.aggregate = function(x) mean(x,na.rm=T))
  
  for(sp in sample_all){
    elementMetadata(motif_in)[[sp]]=as.numeric(NA)
    if(!is.null(stat_in_df_stat[[sp]])){
     
      elementMetadata(motif_in)[[sp]][stat_in_df_stat$qt]=stat_in_df_stat[[sp]]
 
  #add column of NA if motif don't have all samples

    }
    
  }
  print(proc.time()[[3]]-tt1)
  gc()
  return(motif_in)
}
#Find unique type of mutations in SNP analysis in human
unique_mutation<-function(mutation_in){
  mutation_in=data.table(all_SNP=mutation_in)
  mutation_in$unique_class="NA"
  for(i in 1:nrow(mutation_in)){
    SNP_in=mutation_in$all_SNP[i]
    SNP_in_rev=paste0(gsub('.*->','',SNP_in),'->',gsub('->.*','',SNP_in))
    SNP_in2=c(SNP_in,SNP_in_rev)
    if(all(mutation_in[all_SNP%in%SNP_in2]$unique_class=="NA")){
      mutation_in[all_SNP%in%SNP_in2]$unique_class=SNP_in
    }
    
  }
  mutation_out=mutation_in$unique_class
  names(mutation_out)=mutation_in$all_SNP
  return(mutation_out)
}
#Convert trinucleotide SNP into SNP
tri_to_SNP<-function(tri){
  tri_l=gsub('->.*','',tri)
  tri_r=gsub('.*->','',tri)
  return(paste0(substring(tri_l, 2, 2),'->',substring(tri_r, 2, 2)))
}
#Get reverse compliment of SNP or trinucleotide
reverse_comp_SNP<-function(x){
  return(paste0(reverseComplement(DNAString(gsub('->.*','',x))),'->',
                reverseComplement(DNAString(gsub('.*->','',x)))))
  
  
}
#Calculate the dNME based on genotype on left and right of unique genotype and taking into consideration of reverse compliment
dNME_relative_calc<-function(genome1_tri,genome2_tri,NME1,NME2,tri_SNP,tri_SNP_unique,SNP) {
  #Check if trinucleotide is in tri_SNP_unique
  #Getting reverse compliment from original trinucleotide analysis
  genome1_tri_rev=as.character(reverseComplement(DNAString(genome1_tri)))
  genome2_tri_rev=as.character(reverseComplement(DNAString(genome2_tri)))
  #Checking cases: if it's reverse compliment or original
  if(all(c(genome1_tri,genome2_tri) %in% unlist(strsplit(tri_SNP_unique,'->')))){
    return(c(NME1,NME2)[which(c(genome1_tri,genome2_tri)==gsub('->.*','',tri_SNP_unique))]-
             c(NME1,NME2)[which(c(genome1_tri,genome2_tri)==gsub('.*->','',tri_SNP_unique))])
  }  else 
    if(all(c(genome1_tri_rev,genome2_tri_rev)%in%unlist(strsplit(tri_SNP_unique,'->')))){
      return(c(NME1,NME2)[which(c(genome1_tri_rev,genome2_tri_rev)==gsub('->.*','',tri_SNP_unique))]-
               c(NME1,NME2)[which(c(genome1_tri_rev,genome2_tri_rev)==gsub('.*->','',tri_SNP_unique))])
    }else 
      if(all(c(genome1_tri_rev,genome2_tri)%in%unlist(strsplit(tri_SNP_unique,'->')))){
        return(c(NME1,NME2)[which(c(genome1_tri_rev,genome2_tri)==gsub('->.*','',tri_SNP_unique))]-
                 c(NME1,NME2)[which(c(genome1_tri_rev,genome2_tri)==gsub('.*->','',tri_SNP_unique))])
      }else 
        if(all(c(genome1_tri,genome2_tri_rev)%in%unlist(strsplit(tri_SNP_unique,'->')))){
          return(c(NME1,NME2)[which(c(genome1_tri,genome2_tri_rev)==gsub('->.*','',tri_SNP_unique))]-
                   c(NME1,NME2)[which(c(genome1_tri,genome2_tri_rev)==gsub('.*->','',tri_SNP_unique))])
        }
  
}

read.agnostic.mouse<-function(fn,in_dir,replicate="all"){
  fn_sub=gsub('mm10_|.bedGraph|_all_allele_agnostic','',fn)

  #tissue,stage,stat_type,replicate
  tissue=sub('_.*','',fn_sub)
  stat_type=sub('.*_','',fn_sub)
  stage=gsub(paste0(tissue,"_|_",stat_type),'',fn_sub)
  
  file_in=paste0(in_dir,fn)

  cat('processing:',file_in,'\n')
  informME_in=import.bedGraph(file_in)
  if(length(informME_in)>0){
    colnames(elementMetadata(informME_in))=c('score','N','K')
    if(all(seqlevels(informME_in)==gsub('chr','',seqlevels(informME_in)))){seqlevels(informME_in)=paste('chr',seqlevels(informME_in),sep='')}
    #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(informME_in)=start(informME_in)-1
    informME_in$tissue=tissue
    stage=gsub('_5','.5',stage)
    stage=gsub('day','E',stage)
    stage=gsub('E0','P0',stage)
    informME_in$stage=stage
    informME_in$bioreplicate=replicate
    informME_in$Sample=paste(tissue,stage,replicate,sep='-')
    informME_in$stat_type=stat_type
    return(informME_in)
  }
}
#UC
read.agnostic.mouse.uc<-function(file_in,matrix=FALSE,fileter_N=1,gff_in=NA){

  cat('processing:',file_in,'\n')
  informME_in=import.bedGraph(file_in)
  if(length(informME_in)>0){
    colnames(elementMetadata(informME_in))=c('score','N','K')
    if(all(seqlevels(informME_in)==gsub('chr','',seqlevels(informME_in)))){seqlevels(informME_in)=paste('chr',seqlevels(informME_in),sep='')}
    #fit  bedGraph reads, import.bedGraph will remove 1 from start
    start(informME_in)=start(informME_in)-1
    #process file name
    file_in=strsplit(file_in,'\\/')[[1]]
    file_in=file_in[length(file_in)]
    comp= strsplit(strsplit(file_in,'_uc.bedGraph')[[1]],'-vs-')[[1]]
    strain=unlist(lapply(strsplit(comp,'_'),function(x) x[1]))
    #if  contain BL6DBA, use ref is BL6DBA
    strain=ifelse('BL6DBA'%in%strain,'BL6DBA','mm10')
    comp=unlist(lapply(strsplit(comp,'_'),function(x) paste(x[-1],collapse = '_')))
    comp=comp[comp!='']
    comp_stage=unlist(lapply(comp,function(x) {x_split=strsplit(x,'_')[[1]]
    x_split=x_split[-length(x_split)][-1]
    x_split=paste(x_split,collapse = '_')
    return(x_split)}))
    tissue1=strsplit(comp[1],'_')[[1]][1]
    tissue2=strsplit(comp[2],'_')[[1]][1]
    #if BL6DBA, the 1st comp_stage is empty

    comp_stage=gsub('_5','.5',comp_stage)
    comp_stage=gsub('day','E',comp_stage)
    comp_stage=gsub('E0','P0',comp_stage)
    replicate=strsplit(comp[1],'_')[[1]][length(strsplit(comp[1],'_')[[1]])]
    replicate=gsub('merged','',replicate)
    informME_in$Sample=paste0(tissue1,'-',comp_stage[1],'-',tissue2,'-',comp_stage[2],'-',replicate)
    informME_in=informME_in[informME_in$N>=fileter_N]
    informME_in$tissue=tissue1
    informME_in$stage=paste0(comp_stage[1],'-',comp_stage[2])
    informME_in$replicate=replicate
    cat('Minimum N:',min(informME_in$N),'\n')
    #informME_in$Ref=strain
    if(matrix){
      informME_in_dt=as.data.table(mcols(informME_in))[,c("score","Sample")]
      colnames(informME_in_dt)=c("UC","Sample")
      informME_in_dt$UC=as.numeric(informME_in_dt$UC)
      informME_in_dt$region=paste0(seqnames(informME_in),":",start(informME_in),"-",end(informME_in))
      informME_in_dt=informME_in_dt[match(gff_in,region),"UC"]
      colnames(informME_in_dt)=paste0(tissue1,'_',comp_stage[1],'-',tissue2,'_',comp_stage[2],'-',replicate)
      return(informME_in_dt)
    }
    else{return(informME_in)}
  }
}
agnostic_matrix_conversion<-function(gr_in,stat='NME'){
  gr_out=granges(unique(gr_in))
  olap=findOverlaps(gr_in,gr_out,type='equal')
  stat_in_df=elementMetadata(gr_in[queryHits(olap)])[c(stat,'Sample')]
  stat_in_df$idx=NA
  stat_in_df$idx[queryHits(olap)]=subjectHits(olap)
  stat_in_df=as.data.table(stat_in_df)
  
  stat_in_df_stat=dcast.data.table(data=stat_in_df,formula=idx~Sample,value.var = stat,fun.aggregate=mean)#remove agg.fun for new run
  gr_out=gr_out[stat_in_df_stat$idx]
  mcols(gr_out)=stat_in_df_stat[,-1]
  return(gr_out)
  
}
#filtering UC region based on available data
UC_filtering<-function(d){
  d <- sapply(d,function(am) {
    am <- am[,!grepl('P0',colnames(am))]
    am <- am[complete.cases(am),]
  })
  k <- table(unlist(sapply(d,rownames)))
  id <- names(k)[k==length(d)]
  d <- sapply(d,function(i) i[id,],simplify = F)
  return(d)
}

cor_dMML_dNME_enrich<-function(all_region_ts,prob_cutoff,FeDMR){
  all_region_ts$change=FALSE
  cutoff_cor=quantile(all_region_ts$cor,prob=prob_cutoff)
  all_region_ts$change[all_region_ts$cor>=cutoff_cor]=TRUE
  #Find regions that have FeDMR in at least 1 tissue
  
  olap=findOverlaps(all_region_ts,FeDMR)
  all_region_ts$FeDMR=FALSE
  all_region_ts$FeDMR[queryHits(olap)]=TRUE
  cont_table=matrix(c(sum(all_region_ts$change&all_region_ts$FeDMR),
                      sum(all_region_ts$change&(!all_region_ts$FeDMR)),
                      sum((!all_region_ts$change)&(all_region_ts$FeDMR)),
                      sum((!all_region_ts$change)&(!all_region_ts$FeDMR))),nrow=2)
  OR=fisher.test(cont_table)
  print(cont_table)
  return(data.table(OR=OR$estimate,pvalue=OR$p.value,cutoff_cor=cutoff_cor,
                    lowerCI=OR$conf.int[1],upperCI=OR$conf.int[2]))
}

corfunc <- function(m1,m2,type='concordant') {
  if (type=='concordant') {
    rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
  } else {
    scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)            
  }
}
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}
convert_GR<-function(x,direction="GR"){
  if(direction=="GR"){
  strand=unlist(lapply(strsplit(x,','),function(x)x[2]))
  strand=ifelse(is.na(strand),"*",strand)
  x=ifelse(is.na(strand),x,sub(paste0(',\\','-'),'',sub(paste0(',\\','+'),'',x)))
  
  gr=GRanges(seqnames=sub(':.*','',x),
             IRanges(start=as.numeric(sub('-.*','',sub('.*:','',x))),
                     end=as.numeric(sub('.*-','',x))),strand=strand)
  return(gr)}else
    if(direction=="DT"){

      x_dt=as.data.table(mcols(x))
      x_dt$region=paste0(seqnames(x),':',start(x),'-',end(x))
      return(x_dt)
    }else
      if(direction=="matrix"){
        x_mt=as.matrix(mcols(x))
        rownames(x_mt)=paste0(seqnames(x),':',start(x),'-',end(x))
        return(x_mt)
      }
  
}
#Get tss
get_mm10_tss<-function(){
  gtf=fread('../downstream/input/mouse_analysis/grcm38.gtf',data.table=F)
  promoter_in=gtf <- gtf[gtf[,3]=='gene',]
  type <- sub('\".*','',sub('.*gene_type \"','',gtf[,9]))
  gtf <- gtf[type=='protein_coding',]
  gn <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
  gr <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand = gtf[,7])
  names(gr) <- gn
  tss <- promoters(gr,upstream=0,downstream=1)
  return(tss)
  
}
#GO annotation
GO_run<-function(gl,back,cluster,ptcount=0,mapping="org.Mm.eg.db"){
  
  geneList <- factor(as.integer(back %in% gl))
  names(geneList) <- back
  suppressMessages({GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},
                                  annot = annFUN.org, mapping = mapping, ID = "Symbol")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")})
  # pval <- score(resultFisher)
  # pval_adj <- p.adjust(pval, method="BH")
  sigres <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(resultFisher@score),orderBy="classicFisher",numChar=1000)
  sigres<-sigres[sigres$Annotated>=10,]
  #sigres$FDR <- pval_adj[sigres$GO.ID]
  sigres$classicFisher[sigres$classicFisher=="< 1e-30"] <- 0
  sigres$FDR <-p.adjust(sigres$classicFisher,method='BH')
  fc <- ((sigres[,"Significant"])/(sum(GOdata@allScores[GOdata@feasible]==1)))/((sigres[,"Annotated"])/(sum(GOdata@feasible)))
  sigres <- data.frame(sigres,FC=fc)
  sigres <- sigres[order(as.numeric(sigres$FDR),-sigres$FC),]
  sigres=as.data.table(sigres)
  siggene_forID=lapply(sigres$GO.ID,function(x,GOdata){
    gene=sigGenes(GOdata)[sigGenes(GOdata)%in%unlist(genesInTerm(GOdata, x))]
    gl_dt=data.table(rank=1:length(gl),gene=gl)
    mt=match(gl_dt$gene,gene)
    mt=mt[!is.na(mt)]
    highest_rank=NA
    return(list(paste(gene[mt],collapse =";"),highest_rank))
    
    
  },GOdata=GOdata)
  siggene=unlist(lapply(siggene_forID,function(x) x[[1]]))
  max_rank=unlist(lapply(siggene_forID,function(x) x[[2]]))
  if(nrow(sigres)>0){ 
    sigres$genes=siggene
    sigres$higest_ranks=max_rank
    
  }
  sigres$cluster=cluster
  sigres$feasible_allscore=sum(GOdata@allScores[GOdata@feasible]==1)
  sigres$feasible=sum(GOdata@feasible)
  #Find probability of P(sig annotaetd =0)
  
  #sigres$sig_all=length(gl)
  #sigres$bg_all=length(bg)
  #GO data store significant genes
  sigres$sig_all=sum(GOdata@allScores[GOdata@feasible]==1)
  #GO data store total annotated genes
  sigres$bg_all=sum(GOdata@feasible)
  sigres$non_anno_non_sig=sigres$bg_all-sigres$Annotated-sigres$sig_all
  #Annotated non-sig = Annotated
  #non_annotated_sig=total sig
  #Annotated_sig=0
  #choose(a+b,a)=1 given a=0
  sigres$p0_choose=choose(sigres$bg_all-sigres$Annotated,sigres$sig_all)/choose(sigres$bg_all,sigres$sig_all)
  sigres$p0=phyper(0,sigres$sig_all,sigres$bg_all-sigres$sig_all,sigres$Annotated)
  sigres=sigres[Significant!=0]
  sigres$p_fs=sigres[,list(p_fisher=fisher_GO(Significant,Annotated,bg_all,sig_all)),by = seq_len(nrow(sigres))]$p_fisher
  sigres$p_cond=as.numeric(sigres$p_fs)/(1-sigres$p0)
  
  sigres$FDR=p.adjust(sigres$p_cond,method="BH")
  return(sigres)
}
fisher_GO<-function(Significant,Annotated,bg_all,sig_all){
  cont_table=matrix(c(Significant ,sig_all-Significant,Annotated-Significant,bg_all-Annotated-sig_all+Significant),nrow=2)
  return(fisher.test(cont_table)$p.value)
  
  
}
GO_run_tissue<-function(ts,dir_in,enc_type,region_type_sel=NA,bg=NULL,
                         active_enc=F,enc_cor=NA){
  #ranking_stat = "dNME_maxUC" or "dMML_maxUC"
  GO_out_all=list()
  cat("Processing:",ts,'\n')
  fn=paste0(ts,'.csv')
  #read in csv file for given tissue
  csv_in_ts=fread(paste0(dir_in,fn))

  #Note some times Jason use dNME_maxJSD_rank
  csv_in_ts=csv_in_ts[order(dNME_max_UC_pair_adj,decreasing = T)]

  # Getting enhancer
  print(enc_type)
  if(enc_type=="enhancer"&(!active_enc)){
    enhancer=readRDS(bin_enhancer_rds)
    csv_in_gr=convert_GR(csv_in_ts$regions)
   
    mcols(csv_in_gr)=csv_in_ts
    olap=findOverlaps(csv_in_gr,enhancer)
    csv_in_gr=csv_in_gr[queryHits(olap)]
    csv_in_gr$gene=enhancer$`Target Gene`[subjectHits(olap)]
    csv_in_gr$distance=NA
    csv_in_ts=as.data.table(mcols(csv_in_gr))

  }else 
    if(enc_type=="enhancer"&active_enc){
      cat("Analyzing active enhancer\n")
      enc_cor_ts=enc_cor[tissue==ts&cor_FDR<=0.1]
      olap=findOverlaps(convert_GR(csv_in_ts$regions,direction='GR'),
                      convert_GR(enc_cor_ts$region,direction='GR'))
      csv_in_ts=csv_in_ts[queryHits(olap)]
      csv_in_ts$gene=enc_cor_ts[subjectHits(olap)]$target_gene
  }else
    if(enc_type=="promoter"){
      csv_in_ts=csv_in_ts[abs(distance)<=2000]
      
      
    }

  if(region_type_sel!="all"){
    csv_in_ts=csv_in_ts[region_type==region_type_sel]
    print(csv_in_ts)
  }
  #GO annotation
  if(nrow(csv_in_ts)>1){
    #GO annotation for each cluster
    print(nrow(csv_in_ts))
  
    csv_out=lapply(1:max(csv_in_ts$cluster),function(clu){
      sp=paste0(ts,'-',clu)
      csv_in_ts_clu=csv_in_ts[cluster==clu]
      if(nrow(csv_in_ts_clu)>1){
        cat('start processing cluster:',clu,'\n')
        tt1=proc.time()[[3]]
        cat('length of background gene:',length(bg),'\n')
        GO_out_cluster=GO_run(unique(csv_in_ts_clu$gene),bg,cluster=clu)
        csv_in_ts_clu$GO_result=unlist(lapply(csv_in_ts_clu$gene,function(x) paste(GO_out_cluster$Term[grepl(x,GO_out_cluster$genes)],collapse = ';')))
        
        cat('Finish processing cluster:',clu,'in:',proc.time()[[3]]-tt1,'\n')
        return(list(GO_out_cluster_all=GO_out_cluster,csv_in_ts_clu=csv_in_ts_clu))
        
        
        
      }
    })
  }
  print(csv_out[[1]]$GO_out_cluster_all)
  return(csv_out)
}


dcast_matrix<-function(dt_in,value_in,order=T){
  dt_in=dcast.data.table(dt_in,Term~tissue_clu,value.var  = value_in,fill=1)

  dt_in_mt=as.matrix(dt_in[,-1])
  rownames(dt_in_mt)=dt_in$Term
  if(nrow(dt_in_mt)>1&order==T){
  dt_in_mt=dt_in_mt[order(max.col(dt_in_mt)),]
  #dt_in_mt=dt_in_mt[,order(as.numeric(colnames(dt_in_mt)))]
  }
  
  
  return(dt_in_mt)
}
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}
select_top_GO<-function(GO_in,tissue_all,ptcount=0,FDR_cutoff=0.1,FC_cutoff=1.5){
  tissue_all_merged=do.call(rbind,lapply(tissue_all,function (x) {
    
    GO_in_ts=GO_in[[x]]
    
    #Merge all tissue together and recalculate FC with pt count if necessary
    GO_in_ts=do.call(rbind,lapply(GO_in_ts,function(clu_ts){
      
      if(!is.null(clu_ts$GO_out_cluster_all)){
        clu_ts=clu_ts$GO_out_cluster_all
        clu_ts$FC_raw=clu_ts$FC
        clu_ts$FC=((clu_ts$Significant+ptcount)/(clu_ts$feasible_allscore+ptcount))/((clu_ts$Annotated+ptcount)/(clu_ts$feasible+ptcount))
        
        return(clu_ts)
      }
    }))
    GO_in_ts$tissue=x
    return(GO_in_ts)
  }))
  print(head(tissue_all_merged))
  tissue_all_merged$tissue_clu=paste0(tissue_all_merged$tissue,'-',tissue_all_merged$cluster)
  
  #selected top 5 GO terms for each cluster & tissue
  tissue_all_merged_top=do.call(rbind,lapply(unique(tissue_all_merged$tissue_clu),function (x) {
    
    clu_ts=tissue_all_merged[tissue_clu==x]
    #Selection criteria, passing FDR rank by FC, use FDR and p to break ties
    clu_ts_sel=clu_ts[FC>=FC_cutoff&FDR<=FDR_cutoff][order(-FC,FDR,p_cond,decreasing=F)]
    #For cluster with less than 5 significant terms, rank by FC to fill the rest
    if(nrow(clu_ts_sel)<5){
      
      clu_ts_sel=rbind(clu_ts_sel,
                       clu_ts[FC>=FC_cutoff][order(-FC,FDR,p_cond,decreasing=F)][1:(5-nrow(clu_ts_sel))])
    }
    else{clu_ts_sel=clu_ts_sel[1:5]}
    return(clu_ts_sel)
  }))
  return(list(tissue_all_merged=tissue_all_merged,tissue_all_merged_top=tissue_all_merged_top))
  
}
plot_GO_heatmap_all<-function(tissue_all,GO_in,region_type,ptcount=0,FDR_cutoff=0.1,FC_cutoff=1.5,enc_type="enhancer",
                              dir_plot='../downstream/output/mouse_analysis/GO_analysis/kmeans_N17_10run/',plot_pdf=T){
  
  select_top_GO_out=select_top_GO(GO_in,tissue_all,ptcount=ptcount,FDR_cutoff=FDR_cutoff,FC_cutoff=FC_cutoff)
  tissue_all_merged=select_top_GO_out$tissue_all_merged
  tissue_all_merged_top=select_top_GO_out$tissue_all_merged_top
  #Plot for all samples
  figure_dir_all=paste0(dir_plot,'all_sample/')
  if(!dir.exists(figure_dir_all)){dir.create(figure_dir_all)}
  figure_dir_single=paste0(dir_plot,'single_sample/')
  if(!dir.exists(figure_dir_single)){dir.create(figure_dir_single)}
  sel_term=plot_GO_heatmap(
    tissue_all_merged=tissue_all_merged,
    tissue_all_merged_top=tissue_all_merged_top,
    FDR_cutoff=FDR_cutoff,
    fn=paste0(figure_dir_all,'GO_all_samples_',region_type,'_',enc_type,'.pdf')
  )
  #Plot for single sample
  for(ts in unique(tissue_all_merged$tissue)){
    cat("Plotting",ts,'\n')
    tissue_all_merged_top_ts=tissue_all_merged_top[tissue==ts]
    
    tissue_all_merged_top_ts=tissue_all_merged_top_ts[Term %in% sel_term]
    
    
    sel_term_ts=plot_GO_heatmap(
      tissue_all_merged=tissue_all_merged[tissue==ts],
      tissue_all_merged_top=tissue_all_merged_top_ts,
      FDR_cutoff=FDR_cutoff,
      fn=paste0(figure_dir_single,'GO_single_sample_',ts,'_',region_type,'_',enc_type,'.pdf'),
      term_ft=F,plot_pdf=T
    )
    
    
  }  
}


plot_GO_heatmap<-function(tissue_all_merged,tissue_all_merged_top,fn,FDR_cutoff,term_ft=T,plot_pdf=T){
  tissue_all_merged_sel=tissue_all_merged[Term %in% tissue_all_merged_top$Term]
  if(nrow(tissue_all_merged_sel)>0){
    #reshape main matrix and FDR matrix
    # clu_ts_names=unique(tissue_all_merged_sel$tissue_clu)
    clu_ts_names= paste0(expand.grid(1:10,unique(tissue_all_merged$tissue))$Var2,'-',expand.grid(1:10,unique(tissue_all_merged$tissue))$Var1)
    tissue_all_merged_sel_mt=dcast_matrix(tissue_all_merged_sel,"FC",order=T)
    tissue_to_fill=clu_ts_names[!(clu_ts_names %in% colnames(tissue_all_merged_sel_mt))]
    tissue_to_fill_mt=matrix(1,nrow=nrow(tissue_all_merged_sel_mt),ncol=length(tissue_to_fill))
    colnames(tissue_to_fill_mt)=tissue_to_fill
    tissue_all_merged_sel_mt=cbind(tissue_all_merged_sel_mt,tissue_to_fill_mt)
    #There are cases that some cluster have no annotation of any terms due to too few genes in dMML and dNME, we'll fill columns with 1 FC and 1 FDR
    tissue_all_merged_sel_FDR_num=dcast_matrix(tissue_all_merged_sel,"FDR",order=F)
    tissue_all_merged_sel_FDR_num=cbind(tissue_all_merged_sel_FDR_num,tissue_to_fill_mt)
    if(term_ft){
      #Select terms with at least one significant
      sel_term=rownames(tissue_all_merged_sel_mt)
      sel_term=tissue_all_merged_top[
        Term%in%rownames(tissue_all_merged_sel_FDR_num)[which(apply(tissue_all_merged_sel_FDR_num,1,function(x) !all(x>FDR_cutoff)))]]$Term
    }else{
      sel_term=tissue_all_merged_top$Term
    }
    tissue_all_merged_sel_FDR=matrix("",nrow=nrow(tissue_all_merged_sel_FDR_num),ncol=ncol(tissue_all_merged_sel_FDR_num))
    tissue_all_merged_sel_FDR[tissue_all_merged_sel_FDR_num<=FDR_cutoff]="*"
    dimnames(tissue_all_merged_sel_FDR)=dimnames(tissue_all_merged_sel_FDR_num)
    #Add clolumn annotation
    colann= data.frame(cluster=gsub('.*-','',colnames(tissue_all_merged_sel_mt)),tissue=gsub('-.*','',colnames(tissue_all_merged_sel_mt)))
    rownames(colann)=colnames(tissue_all_merged_sel_mt)
    tissue_col=mouse_color()
    cluster_col <- brewer.pal(10,'Set3')
    names(cluster_col)=1:10
    if(ncol(tissue_all_merged_sel_mt)>10){
      gaps_col=seq(10,60,by=10)
    }else(
      gaps_col=NULL
    )
    
    
    if(nrow(tissue_all_merged_sel_FDR)>1){
      tissue_all_merged_sel_mt=scalematrix(tissue_all_merged_sel_mt[sel_term,clu_ts_names])
      tissue_all_merged_sel_FDR=tissue_all_merged_sel_FDR[sel_term,clu_ts_names]
    }
    #It's not useful to plot heatmap with only one column
    if(ncol(tissue_all_merged_sel_mt)>1){
      pheatmap(tissue_all_merged_sel_mt,cluster_rows =F,cluster_cols = F,
               show_colnames = F,show_rownames = T,display_numbers=tissue_all_merged_sel_FDR,border_color = NA,
               color = colorRampPalette(brewer.pal(n = 7, name ="GnBu"))(100),number_color = "white",
               filename=fn,
               cellwidth=60,cellheight=25,annotation_legend = F,angle_col = "90",
               fontsize=30,legend = F,annotation_col = colann,gaps_col = gaps_col,
               annotation_colors = list(tissue=tissue_col,cluster=cluster_col))
    }
    return(sel_term)
  }
}

matrix_conv<-function(dt_in,value.var){
  out_dc=dcast.data.table(dt_in,region~stage,value.var=value.var)
  rn=out_dc$region
  out_dc=as.matrix(out_dc[,-1])
  rownames(out_dc)=rn
  return(out_dc)
}
#Clustering assignment
cluster_assignment<-function(dir_in,dir_out,cluster_region_out_fn,figure_name,UC_merge,
                            cutoffs=0.1,figure_width=2000,figure_height=2000,res=200){
  ifelse(!dir.exists(file.path(dir_out)), dir.create(file.path(dir_out)), FALSE)
  cat('reading in clustering result\n')
  total_run=10
  #Convert into df with major
  cluster_out=list()
  for(fn in c(paste0('uc_',cutoffs,'_',1:10,'.rds'))){
    cluster_in=readRDS(paste0(dir_in,fn))
    for(ts in names(cluster_in)){
      if(fn==paste0('uc_',cutoffs,'_',1,'.rds')){
        cluster_out[[ts]]=data.table(regions=names(cluster_in[[ts]]),cluster_1=cluster_in[[ts]])
      }else{
        cluster_out[[ts]][[paste0("cluster_",gsub(paste0('uc_|.rds|',cutoffs,'_'),'',fn))]]=cluster_in[[ts]][ cluster_out[[ts]]$region]
        
        
      }
    }
    
    
  }
  cat('finding major cluster\n')
  #Find major cluster use cluster_1 as reference
  cluster_out=lapply(cluster_out,function(x){
    for(i in 1:10){
      
      x[[paste0("major_cluster_",i)]]=as.numeric(NA)
      x[[paste0("major_cluster_in_",i)]]=as.numeric(NA)
      #For each cluster, find major cluster
      for(j in 1:10){
        x[cluster_1==j][[paste0("major_cluster_",i)]]= as.numeric(names(which.max(table(x[cluster_1==j][[paste0("cluster_",i)]]))))  
        x[cluster_1==j][[paste0("major_cluster_in_",i)]]=x[cluster_1==j][[paste0("major_cluster_",i)]]==x[cluster_1==j][[paste0("cluster_",i)]]#If in major cluster
      }
    }
    x$percent_cluster_in=rowSums(x[,grepl("major_cluster_in",colnames(x)),with=FALSE])/(total_run)
    return(x)
    
  })
  
  
  # cat('plotting proportion major cluster\n')
  # 
  # pdf(paste0(figure_path,'proportion_run_kmeans_10_all_regions_mm10.pdf'),width=3,height=3)
  # for(ts in names(cluster_out)){
  #   hist(cluster_out[[ts]]$percent_cluster_in,xlab="Proportion of runs in major cluster",main=ts)
  #   
  #   
  # }
  # dev.off()
  # 
  
  # Find regions belong to major cluster ------------------------------------
  cat('assigning minor cluster based on correlation\n')

  cluster_region_out=list()
  for(ts in names(cluster_out)){
    cluster_out_ts=cluster_out[[ts]]
    UC_ts=UC_merge[[ts]][,grepl('UC-',colnames(UC_merge[[ts]]))]

    #Define core clusters
    core_cluster=cluster_out_ts[percent_cluster_in==1]
    core_cluster=core_cluster[,list(regions,cluster_1,percent_cluster_in)]
    core_cluster=cbind(core_cluster,UC_ts[core_cluster$regions,])
    cols=colnames(core_cluster)[grepl(".5-E",colnames(core_cluster))]
    #find patterns of core clusters
    core_cluster_pattern=core_cluster[,lapply(.SD,mean),.SDcols=cols,by=cluster_1]
    core_cluster_pattern=core_cluster_pattern[order(cluster_1)]
    #Find cluster to assign
    cluster_to_assign=cluster_out_ts[percent_cluster_in>=0.5&percent_cluster_in<1]
    
    cluster_to_assign_UC=UC_ts[cluster_to_assign$regions,]
    #Each row is a region, each column is a cluster
    core_cluster_pattern_mt=as.matrix(core_cluster_pattern[,-1])
    rownames(core_cluster_pattern_mt)=core_cluster_pattern$cluster_1
    cor_cluster_out=cor(t(cluster_to_assign_UC),t(core_cluster_pattern_mt))
    #prepare to assign,make sure rows are consistent
    cor_cluster_out=cor_cluster_out[cluster_to_assign$regions,]
    cluster_to_assign$correlation=rowMax(cor_cluster_out)
    cluster_to_assign$cluster=colnames(cor_cluster_out)[apply(cor_cluster_out,1,which.max)]
    #Summary
    cluster_to_assign=cluster_to_assign[,list(regions,cluster,correlation)]
    cluster_to_assign$region_type="noncore_cluster"
    core_cluster$cluster=core_cluster$cluster_1
    core_cluster$correlation=1
    core_cluster$region_type="core_cluster"
    region_out=rbind(core_cluster[,list(regions,cluster,correlation,region_type)],cluster_to_assign[,list(regions,cluster,correlation,region_type)])
    region_out$tissue=ts
    region_out=region_out
    UC_max_ts=UC_merge[[ts]][,grepl('max',colnames(UC_merge[[ts]]))]
    
    UC_max_ts$UC_max_time_adj =gsub(paste0(ts,'-|-all'),'',UC_max_ts$UC_max_time_adj )
    UC_max_ts$UC_max_time  =gsub(paste0(ts,'-|-all'),'',UC_max_ts$UC_max_time  )
    region_out=cbind(region_out,UC_max_ts[region_out$regions,])
    cluster_region_out[[ts]]=region_out
    
    cat("Percent left for:",ts,nrow(region_out)/nrow(cluster_out_ts),'\n')
    write.csv(region_out,paste0(dir_out,ts,'.csv'))
  }
  #UC=0.1 result
  # Percent left for: EFP 0.9988164
  # Percent left for: forebrain 0.9481921
  # Percent left for: heart 0.9686033
  # Percent left for: hindbrain 0.9411712
  # Percent left for: limb 0.9991228
  # Percent left for: liver 0.8832303
  # Percent left for: midbrain 0.9891808
  
  saveRDS(cluster_region_out,cluster_region_out_fn)
  plot_heatmap_cluster(UC_merge,clu=cluster_region_out,figure_name,
  figure_width=figure_width,figure_height=figure_height,res=res)
}
  # Plot heatmap ------------------------------------------------------------
plot_heatmap_cluster<-function(d,clu,figure_name,figure_width=2000,figure_height=2000,res=200){
    cat('Plotting heatmap\n')
    library(gplots)

    d=lapply(d,function(x) x[,!grepl('max',colnames(x))])
    d <- lapply(d,function(x) x[,grepl('UC-',colnames(x))])
    #names(d)=names(UC_merge)
    # dmml <-readRDS(dmml_cor_file)
    # dnme <-readRDS(dnme_cor_file)
    tissue_all=c("EFP","forebrain","heart","hindbrain", "limb","liver" ,"midbrain" )
    timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))


    d=d[tissue_all]
    d=lapply(d,function(x) {
      
      colnames(x)=gsub(paste0('UC-|-all|',paste(tissue_all,'-',sep='',collapse = '|')),'',colnames(x))
      return(x)
    })
    d <- sapply(d,function(i) {
      
      i <- i[rowSums(i) > 0,]
      i <- i[,colnames(i) %in% timeorder]
      i <- i[,order(match(colnames(i),timeorder))]
      
      #i <- scalematrix(i)
      i <- i[complete.cases(i),]
    })
    clu=lapply(clu,function(x){
      out=as.numeric(x$cluster)
      names(out)=x$regions
      return(out)
      
    })
    mat_out=matrix(ncol=39,nrow=0)
    rowann_out=data.frame()
    row_gap=c(0)
    for (n in names(d)) {
      cl <- clu[[n]]
      cl <- sort(cl)
      mat <- do.call(cbind,sapply(tissue_all,function(i) {
        tmp <- matrix(NA,nrow=length(cl),ncol=ncol(d[[i]]),dimnames = list(names(cl),colnames(d[[i]])))
        
        rn <- intersect(names(cl),rownames(d[[i]]))
        tmp[rn,] <- as.matrix(d[[i]][rn,])
        
        colnames(tmp) <- paste0(i,':',colnames(tmp))
        
        tmp
      }))
      na_ma=-which(rowSums(is.na(mat))>0)
      if(length(na_ma)>0){
        mat= mat[na_ma,]
        rowann <- data.frame(tissue_r=n,cluster=sub(':.*','',cl),
                            #dMMLJSDcor=dmml[[n]][rownames(mat)],
                            #dNMEJSDcor=dnme[[n]][rownames(mat)],
                            stringsAsFactors = F)
        rowann=rowann[na_ma,]
      }else{
        
        rowann <- data.frame(tissue_r=n,cluster=sub(':.*','',cl),
                            #dMMLJSDcor=dmml[[n]][rownames(mat)],
                            #dNMEJSDcor=dnme[[n]][rownames(mat)],
                            stringsAsFactors = F)
      }
      #Note for non-tissue specific UC 01, add tissue name to mat
      rownames(mat)=paste0(n,'-',rownames(mat))
      mat_out=rbind(mat_out,mat)
      
      rownames(rowann) <- rownames(mat)
      rowann <- rowann[,ncol(rowann):1]
      rowann_out=rbind(rowann_out,rowann)
      #row_gap=c(row_gap,row_gap[length(row_gap)]+cumsum(rle(sub(':.*','',cl))$lengths))
      row_gap=c(row_gap,row_gap[length(row_gap)]+nrow(mat))
    }
    
    #Refine plotting parameters
    colann <- data.frame(time=sub('.*:','',colnames(mat_out)),tissue=sub(':.*','',colnames(mat_out)),stringsAsFactors = F)
    rownames(colann) <- colnames(mat_out)
    c1 <- mouse_color()
    c2 <- brewer.pal(10,'Set3')
    names(c2) <- 1:10
    c4 <- brewer.pal(length(unique(colann[,1])),'BrBG')
    names(c4) <- sort(unique(colann[,1]))
    #remove row with all NA 

    #sub_sp=sort(sample(1:nrow(mat_out),round(nrow(mat_out))))
      mat_out_sc=scalematrix(mat_out)
      print(head(mat_out_sc))
    png(figure_name,
        width=figure_width,
        height=figure_height,
        res=res,type='cairo')
    if(marker){
        pheatmap(mat_out_sc,cluster_rows = F,
            annotation_row = rowann_out,
            cluster_cols = F,
            annotation_col = colann,show_colnames = F,show_rownames = F,
            #gaps_row = row_gap[-1],
            gaps_col = cumsum(rle(colann[,2])$lengths),
            #filename=figure_name,
            annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4)
                                      #dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10))
                                      
           
    )
    }else{
      pheatmap(mat_out_sc,cluster_rows = F,
            annotation_row = rowann_out,
            cluster_cols = F,
            annotation_col = colann,show_colnames = F,show_rownames = F,
            #gaps_row = row_gap[-1],
            gaps_col = cumsum(rle(colann[,2])$lengths),
            #filename=figure_name,
            annotation_colors = list(tissue=c1,tissue_r=c1,cluster=c2,time=c4)
                                      #dMMLJSDcor=bluered(10),dNMEJSDcor=bluered(10))
                                      
           
    )
    }
    dev.off()
  
}

# correlation analysis ----------------------------------------------------
pval_cor<-function(real_value,perm){
  cor_null=ecdf(perm)
  cor_pval=1-cor_null(real_value)

  return(cor_pval)
}
convert_matrix<-function(perm_in,tissue,stat){
  perm_in_dt=as.data.table(perm_in)
  perm_in_dt$region=rownames(perm_in)
  perm_in_dt$tissue=tissue
  perm_in_dt=melt.data.table(perm_in_dt,id.vars=c('region','tissue'))
  perm_in_dt$variable=paste0(stat,"_perm")
  return(perm_in_dt)
  
}
cor_dt_preprocessing<-function(x,dMML_cor,dNME_cor,dmml_perm,dnme_perm,filtered=FALSE,folder_input=NA) {
  regions=names(dNME_cor[[x]])
  regions=regions[!grepl("X",regions)&!grepl("Y",regions)]
  out_dt=data.table(region=regions,
                    dMML_cor=dMML_cor[[x]][regions],
                    dNME_cor=dNME_cor[[x]][regions],
                    tissue=x)
  
  
  dMML_perm_in=convert_matrix(dmml_perm[[x]],x,'dmml')
  dNME_perm_in=convert_matrix(dnme_perm[[x]],x,'dnme')

  if(filtered==TRUE){  
    csv_in=fread(paste0(folder_input,x,'.csv'))
    out_dt=out_dt[region%in%csv_in$regions]
    out_dt$cluster=csv_in[match(out_dt$region,csv_in$regions)]$cluster
    # out_dt$gene=csv_in[match(out_dt$region,csv_in$regions)]$gene
    # out_dt$distance=csv_in[match(out_dt$region,csv_in$regions)]$distance
    # out_dt$FeDMR=csv_in[match(out_dt$region,csv_in$regions)]$FeDMR
   #Difference between out_dt and csv_in dt from x and y chromosom
    dMML_perm_in=dMML_perm_in[region%in%out_dt$region]
    dNME_perm_in=dNME_perm_in[region%in%out_dt$region]
    enhancer=readRDS(bin_enhancer_rds)
    out_dt_gr=convert_GR(out_dt$region)#20210509, was csv_in$region, bug
    olap=findOverlaps(out_dt_gr,enhancer)
    out_dt$enhancer=FALSE
    out_dt[queryHits(olap)]$enhancer=TRUE
  }
  print(out_dt)
  out_dt$dMML_pval=pval_cor(out_dt$dMML_cor,dMML_perm_in$value)
  out_dt$dNME_pval=pval_cor(out_dt$dNME_cor,dNME_perm_in$value)
  out_dt$dMML_FDR=p.adjust(out_dt$dMML_pval,method='BH')
  out_dt$dNME_FDR=p.adjust(out_dt$dNME_pval,method='BH')
  #dNME_cutoff=min(out_dt[dNME_FDR<=0.25]$dNME_cor)
  #dMML_cutoff=min(out_dt[dMML_FDR<=0.25]$dMML_cor)

  return(out_dt)
  
}
der_calc<-function(x,y){
  der_out=data.table()
  for (i in 1:(length(x)-1)){
    #Drivitative estimation
    der_out=rbind(der_out,
                  data.table(x=x[i],y=y[i],
                             der=(y[i+1]-y[i])/(x[i+1]-x[i])))
    
    
  }
  
  return(der_out)
}
der_flat_finder<-function(der_in,diff_in,density_in,direction=-1,quant=0.05){
  #direction=-1, move to negative, direction=1, move to positive
  der=100
  der_before=100
  #der_max=max(der_in)
  der_quant=quantile(abs(der_in),prob=quant)
  #search from the peak from minimum value, some times a little off from 0
  #(abs(der)!=0)
  #aviod stuck near maximum
  idx=which(density_in==max(density_in[diff_in<=0.05&diff_in>=-0.05]))+20*direction
  while((abs(der)>der_quant |(der_before>der))&(idx %in% 2:(length(der_in)-1))){
    der_before=der_in[idx]
    idx=idx+direction
    #der_p=mean(abs(diff_in[idx:(idx+smooth_window*direction)]/der_max))
    der=der_in[idx]
  }
  return(data.table(der=der,der_before=der_before,der=der_in[idx],idx=idx,x_out=diff_in[idx]))
}
create_folder<-function(folder_out){ifelse(!dir.exists(file.path(folder_out)), dir.create(file.path(folder_out)), FALSE)}
correlation_processing<-function(ts,cor_dt,filtered=F,density_plot=T,FDR_cutoff=0.2,quant=0.25,subsmple_plot=1,
                                 dir_figure,
                                 breaks_color=c("MML only", "NME only", "Both","Neither"),
                                 values_color=c("red", "blue", "green", "purple")){
  theme_density=theme_classic()+theme(legend.position = "bottom",
                                      axis.title.x=element_text(hjust=0.5,size=18,face="bold"),
                                      axis.title.y=element_text(hjust=0.5,size=18,face="bold"),
                                      axis.text.x=element_text(size=16),
                                      axis.text.y=element_text(size=16),
                                      legend.text = element_text(size=16),
                                     )
 
  cat('Processing:',ts,'\n')
  create_folder(dir_figure)
  tissue_in=cor_dt[[ts]]
  #generate MA plot
  tissue_in$cor_diff=tissue_in$dNME_cor -tissue_in$dMML_cor
  tissue_in$cor_mean=(tissue_in$dNME_cor +tissue_in$dMML_cor)/2
  
  #There might be some floating point issue showing mean ==1 but stored as >1
  tissue_in[cor_mean>1]$cor_mean=1
  
  if(density_plot==TRUE){
    cat("Generating density plot\n")
    #pdf(paste0('../downstream/output/correlation/',ts,'_ggplot_raw_density_log_',filtered,'.pdf'),width=4,height=4.5)
     raw_density_log=ggplot(tissue_in,aes(x=dMML_cor,y=dNME_cor))+
      stat_density_2d( geom = "raster",  aes(fill = log(after_stat(density)+0.1)), contour = FALSE,n=200)+
      #scale_fill_distiller(palette = "YlOrBr",trans="log10")
       #scale_fill_distiller(palette = "rainbow",direction = 1)+
       #scale_fill_gradientn(colors=rainbow(20))+
      scale_fill_viridis_c()+theme_density+
     # guides(fill=guide_legend(title="log(density)"))+
      xlab("UC-dMML correlation")+ylab("UC-dNME correlation")
    
    #dev.off()
 
    
   # pdf(paste0('../downstream/output/correlation/',ts,'_raw_log_MA_',filtered,'.pdf'),width=4,height=4.5)
   raw_MA_log=ggplot(tissue_in,aes(x=cor_mean,y=cor_diff))+
      stat_density_2d( geom = "raster",  aes(fill = log(after_stat(density)+0.1)), contour = FALSE,n=200)+
      #scale_fill_distiller(palette = "YlOrBr",trans="log10")
      #scale_fill_distiller(palette = "YlOrRd")+
      scale_fill_viridis_c()+theme_density+
      #guides(fill=guide_legend(title="log(density"))+
      xlab("average correlation")+ylab("correlation difference")
   
    #dev.off()
    #pdf(paste0('../downstream/output/correlation/',ts,'_ggplot_raw_density_',filtered,'.pdf'),width=4,height=4.5)
    raw_density=ggplot(tissue_in,aes(x=dMML_cor,y=dNME_cor))+
            stat_density_2d( geom = "raster",  aes(fill = after_stat(density)), contour = FALSE,n=200)+
            scale_fill_distiller(palette = "YlOrRd",direction = 1)+
            #scale_fill_viridis_c()+theme_density+
            # guides(fill=guide_legend(title="log(density)"))+
            xlab("UC-dMML correlation")+ylab("UC-dNME correlation")
    
    #dev.off()
    
    
    #pdf(paste0('../downstream/output/correlation/',ts,'_raw_MA_',filtered,'.pdf'),width=4,height=4.5)
    raw_MA=ggplot(tissue_in,aes(x=cor_mean,y=cor_diff))+
            stat_density_2d( geom = "raster",  aes(fill = after_stat(density)), contour = FALSE,n=200)+
            #scale_fill_distiller(palette = "YlOrBr",trans="log10")
            scale_fill_viridis_c()+theme_density+
            #guides(fill=guide_legend(title="log(density"))+
            xlab("average correlation")+ylab("correlation difference")
    #dev.off()
 
  }
  
  cat("Finding cutoffs\n")
  tissue_in$cor_mean_round=round(tissue_in$cor_mean*2,digits=1)/2
  cor_cutoffs=data.table()
  pdf(paste0(dir_figure,'cor_cutoff_all_05_',ts,'_all_',filtered,'.pdf'),width=4,height=4.5)
  for(cutoff in sort(unique(tissue_in$cor_mean_round))){
    if(sum(tissue_in$cor_mean_round==cutoff)>=10){
    #Find saddle point for density difference
    tt1=proc.time()[[3]]
    cor_diff_den=density(tissue_in[cor_mean_round==cutoff]$cor_diff,bw="SJ")
    cat('Finish finding density in:',proc.time()[[3]]-tt1,'\n')
    tt1=proc.time()[[3]]
    den_der=der_calc(cor_diff_den$x,cor_diff_den$y)
    cat('Finish derivation calculation:',proc.time()[[3]]-tt1,'\n')
    den_der$diff=den_der$x
    den_der$density=den_der$y
    
    
    #Function to calculate derivatives
    tt1=proc.time()[[3]]
    neg_cutoff=der_flat_finder(den_der$der,den_der$diff,den_der$density,direction=-1,quant=quant)
    pos_cutoff=der_flat_finder(den_der$der,den_der$diff,den_der$density,direction=1,quant=quant)
    cat('Finish finding flat point:',proc.time()[[3]]-tt1,'\n')
    tt1=proc.time()[[3]]
    plot(den_der$diff,den_der$der,xlab="dNME-dMML cor",ylab="1st derivitative",main=cutoff)
    abline(v=pos_cutoff$x_out)
    abline(v=neg_cutoff$x_out)
    plot(cor_diff_den$x,cor_diff_den$y,xlab="dNME-dMML cor",ylab="density",main=cutoff)
    abline(v=pos_cutoff$x_out)
    abline(v=neg_cutoff$x_out)
    cat('Finish ploting derivative:',proc.time()[[3]]-tt1,'\n')
    tt1=proc.time()[[3]]
    #Assigning the catogries based on cutoff
    diff_cutoff_dMML_only=neg_cutoff$x_out
    diff_cutoff_dNME_only=pos_cutoff$x_out
    cutoffs=data.table(x=c(min(tissue_in[cor_mean_round==cutoff]$cor_mean),max(tissue_in[cor_mean_round==cutoff]$cor_mean)),
                       y=c(diff_cutoff_dMML_only,diff_cutoff_dNME_only))
    cat('Finish assigning cutoffs:',proc.time()[[3]]-tt1,'\n')
    tt1=proc.time()[[3]]
    #Plot in density plot
    #Subsample
    subsample_tissue_in=sample(1:nrow(tissue_in),round(nrow(tissue_in)*subsmple_plot),replace = F)
   print(ggplot(tissue_in[subsample_tissue_in],aes(x=cor_mean,y=cor_diff))+
      stat_density_2d( geom = "raster",  aes(fill = after_stat(density)), contour = FALSE,n=200)+
      #scale_fill_distiller(palette = "YlOrBr",trans="log10")
      scale_fill_viridis_c()+theme_density+
      #guides(fill=guide_legend(title="log(density"))+
      xlab("average correlation")+ylab("dNME correlation - dMML correlation")+
      geom_vline(xintercept = cutoffs$x[1],color='red')+
      geom_vline(xintercept = cutoffs$x[2],color='red')+
      geom_line(data=data.table(x=cutoffs$x,y=cutoffs$y[1]),aes(x=x,y=y),color='red')+
      geom_line(data=data.table(x=cutoffs$x,y=cutoffs$y[2]),aes(x=x,y=y),color='red')
   )
    print(ggplot(tissue_in[subsample_tissue_in],aes(x=cor_mean,y=cor_diff))+
      stat_density_2d( geom = "raster",  aes(fill = log(after_stat(density)+0.1)), contour = FALSE,n=200)+
      #scale_fill_distiller(palette = "YlOrBr",trans="log10")
      scale_fill_viridis_c()+theme_density+
      #guides(fill=guide_legend(title="log(density"))+
      xlab("average correlation")+ylab("dNME correlation - dMML correlation")+
      geom_vline(xintercept = cutoffs$x[1],color='red')+
      geom_vline(xintercept = cutoffs$x[2],color='red')+
      geom_line(data=data.table(x=cutoffs$x,y=cutoffs$y[1]),aes(x=x,y=y),color='red')+
      geom_line(data=data.table(x=cutoffs$x,y=cutoffs$y[2]),aes(x=x,y=y),color='red')
    )
    cat('Finish plotting cutoffs:',proc.time()[[3]]-tt1,'\n')
    tt1=proc.time()[[3]]
    cor_cutoffs=rbind(cor_cutoffs,data.table(cor_mean_round=cutoff,
                                             diff_cutoff_dMML_only=diff_cutoff_dMML_only,
                                             diff_cutoff_dNME_only=diff_cutoff_dNME_only))
    
    }
    
  }
  dev.off()
  #Smoothing cutoffs using weighted loess
  cor_cutoffs=cor_cutoffs[order(cor_mean_round)]

  weight_loess= ecdf(tissue_in$cor_mean)(cor_cutoffs$cor_mean_round)
  dMML_fit=loess(diff_cutoff_dMML_only  ~ cor_mean_round,data=cor_cutoffs,span=0.75,
                 weights=weight_loess)
  dNME_fit=loess(diff_cutoff_dNME_only ~ cor_mean_round,data=cor_cutoffs,span=0.75,
                 weights=weight_loess)
  tissue_in$diff_cutoff_dMML_only=predict(dMML_fit,data.table(cor_mean_round=tissue_in$cor_mean))
  tissue_in$diff_cutoff_dNME_only=predict(dNME_fit, data.table(cor_mean_round=tissue_in$cor_mean))
  cat('Finish smooting cutoffs:',proc.time()[[3]]-tt1,'\n')
  tt1=proc.time()[[3]]
  #smoothing plot
  png(paste0(dir_figure,'cor_',ts,'_cluster_all_cor_dMML_sm_quant_weighted_',filtered,'.png'),type='cairo')
  plot(cor_cutoffs$cor_mean_round,cor_cutoffs$diff_cutoff_dMML_only,xlab="mean correlation",ylab="correlation differnce cutoff",xlim=c(-1,1))
  lines(cor_cutoffs$cor_mean_round,predict(dMML_fit,data=cor_cutoffs))
  dev.off()
  png(paste0(dir_figure,'cor_',ts,'_cluster_all_cor_dNME_sm_quant_weighted_',filtered,'.png'),type='cairo')
  plot(cor_cutoffs$cor_mean_round,cor_cutoffs$diff_cutoff_dNME_only,xlab="mean correlation",ylab="correlation differnce cutoff")
  lines(cor_cutoffs$cor_mean_round,predict(dNME_fit,data=cor_cutoffs))
  dev.off()
  cat('Finish plotting smoothed cutoffs:',proc.time()[[3]]-tt1,'\n')
  tt1=proc.time()[[3]]
  #Plot MA and density purly with cutoffs

  raw_MA_log_cutoff_only=raw_MA_log+
    geom_line(data=tissue_in[diff_cutoff_dMML_only<0],aes(x=cor_mean,y=diff_cutoff_dMML_only),color='red',size=0.25)+
    geom_line(data=tissue_in[diff_cutoff_dNME_only>0],aes(x=cor_mean,y=diff_cutoff_dNME_only),color='red',size=0.25)

  #Plot density
 
 
  tissue_in$dMML_cutoff_dMML_only=(tissue_in$cor_mean*2-tissue_in$diff_cutoff_dMML_only)/2
  tissue_in$dNME_cutoff_dMML_only=(tissue_in$cor_mean*2+tissue_in$diff_cutoff_dMML_only)/2
  tissue_in$dMML_cutoff_dNME_only=(tissue_in$cor_mean*2-tissue_in$diff_cutoff_dNME_only)/2
  tissue_in$dNME_cutoff_dNME_only=(tissue_in$cor_mean*2+tissue_in$diff_cutoff_dNME_only)/2
 #Make sure the cutoffs is in range
  raw_density_log_cutoff_only=raw_density_log+
    geom_point(data=tissue_in[diff_cutoff_dMML_only<0&
                                dMML_cutoff_dMML_only >-1&
                                dMML_cutoff_dMML_only <=1],aes(x=dMML_cutoff_dMML_only,y=dNME_cutoff_dMML_only),color='red',size=0.25)+
    geom_point(data=tissue_in[diff_cutoff_dNME_only>0&
                              dMML_cutoff_dNME_only >-1&
                                dMML_cutoff_dNME_only <=1],aes(x=dMML_cutoff_dNME_only,y=dNME_cutoff_dNME_only),color='red',size=0.25)
  #Plot density and FDR with cutoffs
  cutoff_dMML=min(tissue_in[dMML_FDR<=FDR_cutoff]$dMML_cor)
  cutoff_dNME=min(tissue_in[dNME_FDR<=FDR_cutoff]$dNME_cor)
  if(is.infinite(cutoff_dNME)){cutoff_dNME=1}
  if(is.infinite(cutoff_dMML)){cutoff_dMML=1}
  dMML_seq=c(seq(-1,cutoff_dMML,0.001),cutoff_dMML)
  dNME_seq=c(seq(-1,cutoff_dNME,0.001),cutoff_dNME)
  FDR_cutoff_dt=rbind(
    data.table(
      region_type="dMML FDR cutoff",
      dMML=rep(cutoff_dMML,length(dNME_seq)),
      dNME=dNME_seq
    ),
    data.table(
      region_type="dNME FDR cutoff",
      dMML=dMML_seq,
      dNME=rep(cutoff_dNME,length(dMML_seq))
    ))
  raw_density_log_cutoff_FDR=raw_density_log+
    geom_point(data=tissue_in[(dNME_cutoff_dMML_only>cutoff_dNME|dMML_cutoff_dMML_only>cutoff_dMML)&(!(dNME_cutoff_dMML_only>0&dMML_cutoff_dMML_only<0))],
               aes(x=dMML_cutoff_dMML_only,y=dNME_cutoff_dMML_only),color='red',size=0.25)+
    geom_point(data=tissue_in[(dNME_cutoff_dNME_only>cutoff_dNME|dMML_cutoff_dNME_only>cutoff_dMML)&(!(dNME_cutoff_dNME_only>0&dMML_cutoff_dNME_only<0))],
               aes(x=dMML_cutoff_dNME_only,y=dNME_cutoff_dNME_only),color='red',size=0.25)+
    geom_line(data=FDR_cutoff_dt[region_type=='dMML FDR cutoff'],aes(x=dMML,y=dNME),color='red',size=0.5)+
    geom_line(data=FDR_cutoff_dt[region_type=='dNME FDR cutoff'],aes(x=dMML,y=dNME),color='red',size=0.5)
  cat('Finish plotting smoothed cutoffs:',proc.time()[[3]]-tt1,'\n')
  #Assigning regions
  rm(diff_cutoff_dNME_only)
  rm(diff_cutoff_dMML_only)

  tissue_in$region_type="NA"
  #Both
  tissue_in[cor_diff<diff_cutoff_dNME_only&
              cor_diff>diff_cutoff_dMML_only&
              (dNME_FDR<=FDR_cutoff|dMML_FDR<=FDR_cutoff)]$region_type="Both"
  #One significant the other one <0 are also quantified as only
  #dMML_only
  tissue_in[(dNME_FDR<=FDR_cutoff|dMML_FDR<=FDR_cutoff)&
              (cor_diff<=diff_cutoff_dMML_only|
                 (dMML_cor>0&dNME_cor<0))]$region_type="MML only"
  #dNME_only
  tissue_in[(dNME_FDR<=FDR_cutoff|dMML_FDR<=FDR_cutoff)&
              (cor_diff>=diff_cutoff_dNME_only|
                 (dNME_cor>0&dMML_cor<0))]$region_type="NME only"
  
  #Neither
  tissue_in[dNME_FDR>FDR_cutoff&dMML_FDR>FDR_cutoff]$region_type="Neither"
  cat('Finish assign regions types cutoffs:',proc.time()[[3]]-tt1,'\n')
  tt1=proc.time()[[3]]
  #dev.off()
  #Smooth cutoffs
  subsample_tissue_in=sample(1:nrow(tissue_in),round(nrow(tissue_in)*subsmple_plot),replace = F)
  #png(paste0('../downstream/output/correlation/cor_',ts,'_cluster_all_cor_smoothed_quant_weighted_',filtered,'.png'),width=4,height=4.5,units = 'in',res=144)
   catotry_dot_cor=ggplot(tissue_in[subsample_tissue_in],aes(x=dMML_cor,y=dNME_cor,color=region_type,fill=region_type))+
          geom_point(alpha=0.1)+xlab("UC-dMML correlation")+ylab("UC-dNME correlation")+
          guides(colour = guide_legend(override.aes = list(alpha = 1),nrow=2,byrow=TRUE))+
    theme_density+theme(legend.title=element_blank())
  #dev.off()
  #MA plot
  #png(paste0('../downstream/output/correlation/cor_',ts,'_cluster_all_cor_MA_smoothed_quant_',filtered,'.png'),width=4,height=4.5,units = 'in',res=144)
   tissue_in$region_type=factor(tissue_in$region_type, levels=c("MML only","NME only","Both","Neither"))
  catotry_dot_MA=ggplot(tissue_in[subsample_tissue_in],aes(x=cor_mean ,y=cor_diff ,color=region_type,fill=region_type))+
          geom_point(alpha=0.1)+xlab("average correlation")+ylab("correlation difference")+
    guides(colour = guide_legend(override.aes = list(alpha = 1),nrow=2,byrow=TRUE))+
    theme_density+scale_color_manual(breaks = breaks_color, values=values_color)
  #dev.off()
cat('Finish plot dot plot:',proc.time()[[3]]-tt1,'\n')
tt1=proc.time()[[3]]
  
  if(density_plot==TRUE){
      cat("Generating density plot after finish\n")
     #cutoff_dt=dMML_dNME_cutoff_dt(tissue_in,FDR_cutoff)
     #print(head(cutoff_dt))
     mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(8)
    
      png(paste0(dir_figure,ts,'_density_',filtered,'.png'),width=9,height=12.5,units = 'in',res=1080,type='cairo')
      print(ggarrange(raw_density_log,raw_MA_log,raw_density_log_cutoff_only,raw_MA_log_cutoff_only,
                      raw_density_log_cutoff_FDR,
                      #cat_density_log,cat_density_MA_log, 
                      ncol = 2, nrow = 3,common.legend = TRUE))
      dev.off()
      png(paste0(dir_figure,ts,'_dot_',filtered,'.png'),width=4.5,height=9,units = 'in',res=1080,type='cairo')
      print(ggarrange(catotry_dot_cor,catotry_dot_MA, 
                      ncol = 1, nrow = 2,common.legend = TRUE))
      dev.off()
      #cat("Finish generating density plot in:",proc.time()[[3]]-tt1,'\n')
      cat('Finish plot dot plot:',proc.time()[[3]]-tt1,'\n')
      tt1=proc.time()[[3]]
    #}else{cat("high FDR threshold for density plot\n")}
  }
  return(tissue_in)
  
}
dMML_dNME_cutoff_dt<-function(dt_in,FDR_cutoff){
  #Find dMML and dNME cutoff for dMML only and dNME only region: diff=dNME-dMML, dNME=(diff+mean)/2, dMML=(mean-diff)/2
  dt_in$dMML_cutoff_dMML_only=(dt_in$cor_mean*2-dt_in$diff_cutoff_dMML_only)/2
  dt_in$dNME_cutoff_dMML_only=(dt_in$cor_mean*2+dt_in$diff_cutoff_dMML_only)/2
  dt_in$dMML_cutoff_dNME_only=(dt_in$cor_mean*2-dt_in$diff_cutoff_dNME_only)/2
  dt_in$dNME_cutoff_dNME_only=(dt_in$cor_mean*2+dt_in$diff_cutoff_dNME_only)/2
  #reformat into a data.table
  cutoff_dMML=min(dt_in[dMML_FDR<=FDR_cutoff]$dMML_cor)
  dMML_seq=c(seq(-1,cutoff_dMML,0.001),cutoff_dMML)
  cutoff_dNME=min(dt_in[dNME_FDR<=FDR_cutoff]$dNME_cor)
  dNME_seq=c(seq(-1,cutoff_dNME,0.001),cutoff_dNME)
  #FDR cutoffs
  cutoff_dt=rbind(
    data.table(
      region_type="dMML FDR cutoff",
      dMML=rep(cutoff_dMML,length(dNME_seq)),
      dNME=dNME_seq
    ),
    data.table(
      region_type="dNME FDR cutoff",
      dMML=dMML_seq,
      dNME=rep(cutoff_dNME,length(dMML_seq))
    ),
    data.table(
      region_type="dMML only",
      dMML=unique(round(dt_in[region_type=='dMML_only']$dMML_cor,digits=2)),
      dNME=dt_in[region_type=='dMML_only',list(dNME_max=max(dNME_cor)),by=round(dMML_cor,digits=2)]$dNME_max
      
    ),
    data.table(
      region_type="dNME only",
      dMML=dt_in[region_type=='dNME_only',list(dMML_max=max(dMML_cor)),by=round(dNME_cor,digits=2)]$dMML_max,
      dNME=unique(round(dt_in[region_type=='dNME_only']$dNME_cor,digits=2))
      
    )
  )
  cutoff_dt=cutoff_dt[order(dMML)]
  return(cutoff_dt)
}

plot_correlation<-function(tissue_out_filtered,pdf_fn,plot_pdf=T,
                                 breaks_color=c("MML only", "NME only", "Both","Neither"),
                                 values_color=c("red", "blue", "green", "purple")){
  tissue_out_filtered=do.call(rbind,tissue_out_filtered)
  tissue_out_filtered$region_type=factor(tissue_out_filtered$region_type,
                                         levels=c("MML only","NME only","Both","Neither"))
  #All regions
  tissue_out_filtered_frequency=tissue_out_filtered[tissue!="NT",list(count=length(region)),by=list(tissue,region_type)]
  tissue_out_filtered_frequency_p=tissue_out_filtered_frequency[,list(percentage=round(count/sum(count)*100,digits=0),region_type=region_type),by=list(tissue)]
  tissue_out_filtered_frequency_p[,list(minp=min(percentage),maxp=max(percentage)),by=list(region_type)]
  print(t.test(tissue_out_filtered_frequency_p[region_type=="MML only"]$percentage,
         tissue_out_filtered_frequency_p[region_type=="NME only"]$percentage,alternative = "less"))
  plot_out=ggplot(tissue_out_filtered_frequency, aes(y=count, x=tissue,fill=region_type)) + 
          geom_bar( stat="identity",position="fill")+ylab("")+xlab("")+
          
          theme_classic()+theme(axis.text.x=element_text(size=32,angle=90),
                                axis.text.y=element_text(size=32))+
          scale_fill_manual(breaks=breaks_color,values=values_color)+guides(fill=guide_legend(nrow=2,byrow=TRUE))+
          theme(legend.position = "bottom",legend.title = element_blank(),legend.text =element_text(size=28))
  if(plot_pdf){
  pdf(pdf_fn,width=7,height=16)
  print(plot_out)
  dev.off()
  }else{return(plot_out)}
  
}
assign_regions<-function(tissue_out_filtered,folder_in_clu,DNAase){
  #assign regions 

  DNAase=convert_GR(DNAase,direction="DT")
  lapply(tissue_out_filtered,function(cor_dt_in){
    tissue_in=unique(cor_dt_in$tissue)
    csv_in=fread(paste0(folder_in_clu,tissue_in,'.csv'))[,-1]
    print(head(csv_in))
    csv_in$region_type=cor_dt_in[match(csv_in$regions,cor_dt_in$region)]$region_type
    csv_in$DNAase=FALSE
    DNAase_region=which(csv_in$regions%in% DNAase$region)
    csv_in[DNAase_region]$DNAase=TRUE
    tss=get_mm10_tss()
    dt_nearest=GenomicRanges::distanceToNearest(convert_GR(csv_in$regions),tss)
    csv_in$distance=mcols(dt_nearest)$distance
    csv_in$gene=names(tss)[subjectHits(dt_nearest)]
    write.csv(csv_in,paste0(folder_in_clu,tissue_in,'.csv'),row.names = F)
    return(NULL)
  })
  
}
#For promoter enhancer comparison
diff_stat_extraction<-function(UC_merge_in,stat,csv_in){
  
  dt_out=UC_merge_in[rn %in% csv_in$regions,.SD,.SDcols=c(1,which(grepl(paste0(stat,'-'),colnames(UC_merge_in))))]
  dt_out=melt.data.table(dt_out,variable.name='sample',value.name=stat)
  dt_out$tissue=tissue_in
  dt_out$enhancer=csv_in[match(dt_out$rn,regions)]$enhancer
  dt_out$distance=csv_in[match(dt_out$rn,regions)]$distance
  dt_out[,states:="NA"]
  dt_out[enhancer==TRUE,states:="enhancers"]
  dt_out[abs(distance)<=2000,states:="promoters"]
  return(dt_out)
  
}
GO_sheets<-function(GO_result,enc_type,dNME_cor=dNME_cor,dMML_cor=dMML_cor,FeDMR_dir='../downstream/input/FeDMR_Ecker/',
                    motif_Ken_dir='../downstream/input/Ken_motif_binding_site/',FDR_cutoff=0.1,mm10_CpG=cgs,out_dir='../downstream/output/mouse_analysis/GO_analysis/GO_sheets/'){
  #Annotate cluster
  
  #est_time_table=fread('../downstream/input/mouse_analysis/clustering/tissue_specific/heatmap_uc_01/uc_01_tissue_clu_anno.csv')
  
  GO_out=list()
  for(GO_type in names(GO_result)){
    GO_out_result=data.table()
    for(ts in names(GO_result[[GO_type]])){
      cat("Processing:",ts,' ',GO_type,'\n')
      #For each cluster
      gene_sig_ts=do.call(rbind,lapply(GO_result[[GO_type]][[ts]],function(x){
        #Find significant terms
 
        if(!is.null(x$GO_out_cluster_all)){
        GO_sig=x$GO_out_cluster_all[FDR<=FDR_cutoff]
        # if(nrow(GO_sig)>0){
        #   
        #   gene_sig=x$csv_in_ts_clu[gene%in% unique(unlist(strsplit(GO_sig$genes,';'))),
        #                            list(gene,cluster,region,distance,dMML_maxUC,dNME_maxUC,dMML_maxpair,dNME_maxpair)]
        #   # gene_sig$dMML_max_time=UC_merge[[ts]][gene_sig$region,"dMML_max_time"]
        #   # gene_sig$dNME_max_time=UC_merge[[ts]][gene_sig$region,"dNME_max_time"]
        #   # gene_sig$UC_max_time=UC_merge[[ts]][gene_sig$region,"UC_max_time"]
        #   
        #   gene_sig$GO_result=unlist(lapply(gene_sig$gene,function(g) paste(GO_sig[grepl(g,GO_sig$gene)]$Term,collapse=';')))
        #   gene_sig$cluster=unique(x$GO_out_cluster_all$cluster)
        #   return(gene_sig)
        # }
        #return(GO_sig[,list(GO.ID,Term,Annotated,Significant,Expected,FC,p_cond,FDR,genes)])
        return(GO_sig[,list(GO.ID,Term,cluster,FDR,genes)])
        }
      }))
      if(nrow(gene_sig_ts)>0&length(gene_sig_ts)>0&!is.null(gene_sig_ts)){
      write.xlsx(gene_sig_ts,paste0(out_dir,GO_type,'.xlsx'),row.names=F,sheet=ts,append=T)
      }
    }
    #   if(!is.null(gene_sig_ts)&length(gene_sig_ts)>0){
    #   gene_sig_ts$tissue=ts
    #   gene_sig_ts$est_stage=est_time_table[Tissue==ts][match(gene_sig_ts$cluster,Cluster)]$Stage
    #   #Get Ecker's result put in sheets
    #   cat("Assigning FeDMR\n")
    #   feDMR_in_gr_all=GRanges()
    # 
    #     # feDMR_in_all=data.table()
    #     # 
    #     # hg19tomm10=import.chain('../downstream/input/hg19ToMm10.over.chain')
    #     # for(fn in dir(FeDMR_dir,pattern=ts)){
    #     #   feDMR_in=fread(paste0(FeDMR_dir,fn))
    #     #   feDMR_in=feDMR_in[chrom!="chrX"]
    #     #   #double check if this is mm10
    #     #   feDMR_in$region=paste0(feDMR_in$chrom,":",feDMR_in$start,"-",feDMR_in$end)
    #     # 
    #     #   #feDMR_in_gr=unlist(liftOver(feDMR_in_gr,hg19tomm10))
    #     #   stage=gsub('feDMR_','',gsub(paste0('_',ts,'.tsv'),'',fn))
    #     # 
    #     #   feDMR_in$stage=stage
    #     #   feDMR_in_all=rbind(feDMR_in_all,feDMR_in)
    #     # 
    #     # }
    #     # if(!is.null(feDMR_in_all)&nrow(feDMR_in_all)>0){
    #     #   feDMR_in_all=feDMR_in_all[,list(stage=paste(stage,collapse = ';')),by=region]
    #     #   olap=findOverlaps(convert_GR(gene_sig_ts$region),convert_GR(feDMR_in_all$region))
    #     #   gene_sig_ts$FeDMR="NA"
    #     #   gene_sig_ts$FeDMR[queryHits(olap)]=feDMR_in_all[subjectHits(olap)]$stage
    #     # }
    #     cat("Assigning Motif\n")
    #     gene_sig_ts$Ken_dNME=""
    #     gene_sig_ts$Ken_dMML=""
    #     gene_sig_ts$Ken_dNME_CpG=FALSE
    #     motif_Ken_dir='../downstream/input/mouse_analysis/motif_analysis/Ken_motif_locus/'
    #     motif_in_dNME=readRDS(paste0(motif_Ken_dir,ts,"_motif_site_dNME.rds"))
    #     if(length(motif_in_dNME)>0){
    #     motif_in_dNME=Ken_motif_merge(motif_in_dNME)
    #     #This step will only get one motif for each region, however, it's possible there're multiple motifs
    #     gene_sig_ts$Ken_dNME=convert_GR(motif_in_dNME,direction="DT")[match(gene_sig_ts$region,region)]$motif
    #     motif_locus_in=readRDS(paste0(Ken_motif_locus,ts,'_motif_site_dNME_locus.rds'))
    #     motif_locus_in=Ken_motif_merge(motif_locus_in)
    #     print(motif_locus_in)
    #     olap_CG=findOverlaps(motif_locus_in,mm10_CpG,maxgap = 2)
    #     motif_locus_in$CpG_mm10=FALSE
    #     motif_locus_in$CpG_mm10[unique(queryHits(olap_CG))]=TRUE
    #     gene_sig_ts$Ken_dNME_CpG=as.data.table(mcols(motif_locus_in))[match(gene_sig_ts$region,region)]$CpG_mm10
    #     }
    #     motif_in_dMML=readRDS(paste0(motif_Ken_dir,ts,"_motif_site_dMML.rds"))
    #     if(length(motif_in_dMML)>0){
    #     motif_in_dMML=Ken_motif_merge(motif_in_dMML)
    #     gene_sig_ts$Ken_dMML=convert_GR(motif_in_dMML,direction="DT")[match(gene_sig_ts$region,region)]$motif
    #     }
    #     
    #     gene_sig_ts=gene_sig_ts[,list(gene,tissue,region,cluster,est_stage,dMML_maxUC,dNME_maxUC,Ken_dNME,Ken_dMML,Ken_dNME_CpG)]
    #     gene_sig_ts$dNME_cor=dNME_cor[[ts]][gene_sig_ts$region]
    #     gene_sig_ts$dMML_cor=dMML_cor[[ts]][gene_sig_ts$region]
    #     gene_sig_ts=gene_sig_ts[order(dNME_maxUC,round(dNME_cor,digits=1),decreasing = T )]
    #     write.csv(gene_sig_ts,paste0('../downstream/output/GO_sheets/',ts,'_',GO_type,"_",enc_type,'.csv'),row.names=F)
    #     gene_sig_ts=gene_sig_ts[,list(gene,tissue,cluster,est_stage,dMML_maxUC,dNME_maxUC,dNME_cor,dMML_cor,Ken_dNME,Ken_dMML,Ken_dNME_CpG)]
    #     GO_out_result=rbind(GO_out_result,gene_sig_ts)
    #   }
    #  
    # }
    
    # write.csv(GO_out_result,paste0('../downstream/output/mouse_analysis/GO_analysis/GO_sheets/summary sheet/',GO_type,"_",enc_type,'.csv'),row.names=F)
  }
  
}
Ken_motif_merge<-function( motif_Ken_raw){
  
  return(do.call(c,lapply(names(motif_Ken_raw),function(x){
    motif_Ken_raw_in=motif_Ken_raw[[x]]
    motif_Ken_raw_in$motif_full=x
    motif_Ken_raw_in$motif=gsub('.*_','',motif_Ken_raw_in$motif_full)
    return(motif_Ken_raw_in)
    
  })))
  
}
motif_sig_Ken<-function(ts,stat,motif_locus_ken,enhancer_regions_ts,dir_in_Ken='../downstream/input/mouse_analysis/motif_analysis/Ken_motif_result_all_regions/'){
  motif_sig=fread(paste0(dir_in_Ken,ts,'_OR_residual_',stat,'.csv'))
  if(nrow(motif_sig)>0){
    motif_locus=do.call(c,lapply(motif_sig$motif,function(x) {
      motif_out=motif_locus_ken[[ts]][[x]]
      
      motif_out$motif=gsub('.*_','',x)
      return(motif_out)
      
    }))
    
    motif_locus_dt=unique(as.data.table(mcols(motif_locus)))
    
    idx=match(motif_locus_dt$region,enhancer_regions_ts$region)
    
    motif_locus_dt=cbind(motif_locus_dt,enhancer_regions_ts[idx,-1,with=T])[!is.na(gene)]
  
    write.csv(motif_locus_dt,
              paste0('../downstream/output/mouse_analysis/motif_analysis/enhancer_gene_motif/motif_gene_',ts,'_',stat,'.csv'),
             row.names = F)
    motif_locus_dt_collapse=motif_locus_dt[,list(motif=paste(motif,collapse=";")),by=region]
    for(mt in unique(motif_locus$motif)){
      bed_fn=paste0('../downstream/output/mouse_analysis/motif_analysis/motif_locus_bed/',stat,'/',ts,'_',mt,'.bed')
      bed_fn=gsub('::','-',bed_fn)

    write.table(as.data.table( motif_locus)[motif==mt,list(seqnames,start,end,motif)],
                bed_fn,
                row.names = F,col.names = F,quote =F)
    }
    motif_locus_dt_collapse$tissue=ts
    
  }else{
    motif_locus_dt_collapse=data.table(region=NA,motif=NA)
    
    motif_locus=NA
  }
  return(list(motif_locus_dt=motif_locus_dt_collapse,motif_locus=motif_locus))
}
pubmed_rec<-function(genes,keys){
 
  my_query=paste0('(',genes,'[Title/Abstract]) AND ',keys,' AND (Human OR Mouse)')
  print(my_query)
   tt1=proc.time()[[3]]
  my_entrez_id = get_pubmed_ids(my_query)
  
  if(as.numeric(my_entrez_id $Count)>0){
    
    my_abstracts_xml <-do.call(rbind,lapply(articles_to_list(fetch_pubmed_data(pubmed_id_list = my_entrez_id)),function(x) {
     
      PMID= unlist(custom_grep(x,tag="PMID"))
      if(is.null(PMID)){PMID=""}
      title=unlist(custom_grep(x,tag="ArticleTitle"))
      if(is.null(title)){title=""}
      return(data.table(PMID= PMID,
                        #journal_abb=custom_grep(x,tag="ISOAbbreviation"),
                        title=title
                        #year=paste(unique(custom_grep(x,tag="Year")),collapse=";")
                        
      ))
      
      
    }))
    if(!is.null(my_abstracts_xml)){
    
    my_abstracts_xml$gene=genes
    cat("Finish processing", genes, "in ",proc.time()[[3]]-tt1,'\n')
   
    return(my_abstracts_xml)
    }
  }
  
}
OR_repeats<-function(repeat_type_in,re_web_class_type,ts_region,non_ts_region,tissue){
  repeat_class=re_web_class_type[repeat_type==repeat_type_in]$repeat_class
  #Do not use grepl, some name contain other names, like L1MdFanc_I and L1MdFanc_II
  ts_region_repeat=sum(do.call(c,lapply(strsplit(ts_region$repeat_type,';'),function(x) repeat_type_in %in% x)))
  ts_region_non_repeat=length(ts_region)-ts_region_repeat
  non_ts_region_repeat=sum(do.call(c,lapply(strsplit(non_ts_region$repeat_type,';'),function(x) repeat_type_in %in% x)))
  non_ts_region_non_repeat=length(non_ts_region)-non_ts_region_repeat
  OR=fisher.test(matrix(c(ts_region_repeat,ts_region_non_repeat,non_ts_region_repeat,non_ts_region_non_repeat),nrow=2))
  return(data.table(
    OR=OR$estimate,
    pvalue=OR$p.value,
    tissue=tissue,
    repeat_type=repeat_type_in,
    repeat_class=repeat_class,
    observed=ts_region_repeat,
    ts_region_non_repeat=ts_region_non_repeat,
    non_ts_region_repeat,
    non_ts_region_non_repeat,
    expected=non_ts_region_repeat/length(non_ts_region)*length(ts_region)))
}
repeat_olap_all<-function(ts_region,non_ts_region,re_web){
  ts_region=add_repeats(ts_region,re_web)
  non_ts_region=add_repeats(non_ts_region,re_web)
  ts_repeat=sum(ts_region$repeat_type!="NA")
  ts_non_repeat=sum(ts_region$repeat_type=="NA")
  non_ts_repeat=sum(non_ts_region$repeat_type!="NA")
  non_ts_non_repeat=sum(non_ts_region$repeat_type=="NA")
  return(data.table(
    ts_repeat=ts_repeat,
    ts_non_repeat=ts_non_repeat,
    non_ts_repeat=non_ts_repeat,
    non_ts_non_repeat=non_ts_non_repeat,
    total_ts_reion=length(ts_region)
  ))
  
}
repeat_olap_individual<-function(ts_region,non_ts_region,re_web,re_web_class_type,tissue,ncores=12){
  tt1=proc.time()[[3]]
  ts_region=add_repeats(ts_region,re_web)
  non_ts_region=add_repeats(non_ts_region,re_web)
  cat("Selected UC region overlap with repeats:",sum(ts_region$repeat_type!="NA"),"\n")
  cat("Percent selected UC region overlap with repeats:",sum(ts_region$repeat_type!="NA")/length(ts_region)*100,"%\n")
  repeat_all_ts=unique(unlist(strsplit(ts_region$repeat_type,';')))
  repeat_all_ts=repeat_all_ts[repeat_all_ts!="NA"]
  cat("start processing each motif\n")
  repeats_output_ts=mclapply(repeat_all_ts,OR_repeats,
                             re_web_class_type=re_web_class_type,
                             ts_region=ts_region,
                             non_ts_region=non_ts_region,tissue=tissue,
                             mc.cores=ncores)
  repeats_output_ts=do.call(rbind,repeats_output_ts)
  repeats_output_ts$FDR=p.adjust(repeats_output_ts$pvalue,method="BH")

  gc()
  cat("Finish processing in ",proc.time()[[3]]-tt1,'\n')
  return(repeats_output_ts)
}
add_repeats<-function(region_in,repeats_in){
  #add minOverlap here
  min_olap=0.5*mean(width(region_in))
  cat('Min olap setting:',min_olap,'\n')
  region_in_olap=as.data.table(findOverlaps(region_in,repeats_in,minoverlap = min_olap))
  region_in_olap$repeat_type=repeats_in$repeat_type[region_in_olap$subjectHits]
  region_in_olap$repeat_class=repeats_in$repeat_class[region_in_olap$subjectHits]
  region_in_olap=region_in_olap[,list(repeat_type=paste(repeat_type,collapse = ";"),
                                      repeat_class=paste(repeat_class,collapse = ";")),
                                by=queryHits]
  
  region_in$repeat_type="NA"
  region_in$repeat_class="NA"
  region_in$repeat_type[region_in_olap$queryHits]=region_in_olap$repeat_type
  region_in$repeat_class[region_in_olap$queryHits]=region_in_olap$repeat_class
  return(region_in)
}
getCpgSitesmm10 <- function(chrsOfInterest=paste("chr",1:19,sep="")){
  #1:19 exclude the sex chromosome
  # Obtain all CpG sites
  cgs <- lapply(chrsOfInterest, function(x)  GRanges(x,IRanges(start(matchPattern("CG", Mmusculus[[x]])),with=2)))
  # Set genome and seqlengths
  cgs <- setGenomeLengths(do.call('c',cgs),chrsOfInterest=chrsOfInterest,genome_in="mm10")
  # Return CpG site GR
  return(cgs)
}
getCpgSiteshg19 <- function(chrsOfInterest=paste("chr",1:21,sep="")){
  # Obtain all CpG sites
  cgs <- lapply(chrsOfInterest, function(x)  GRanges(x,IRanges(start(matchPattern("CG", Hsapiens[[x]])),with=2)))
  # Set genome and seqlengths
  cgs <- setGenomeLengths(do.call('c',cgs),chrsOfInterest=chrsOfInterest,genome_in="hg19")
  # Return CpG site GR
  return(cgs)
}

#Function for ChIP overlap looking for examples
ChIP_assignment<-function(region_in,factor_in,stat_in="dNME"){
  region_in=region_in[,list(region,tissue,cluster,gene,region_type,UC_max_time,dNME_max_pair,dMML_max_pair,UC_max_pair,dNME_motif,dMML_motif,gene,PMID)]
  factor_in=makeGRangesFromDataFrame(factor_in,keep.extra.columns = T)
  olap=findOverlaps(convert_GR(region_in$region),factor_in)
  region_in_factor=region_in[queryHits(olap)]
  region_in_factor$motif_ChIP=factor_in$metadata[subjectHits(olap)]
  if(nrow(region_in_factor)>0){
  region_in_factor$predict_in_ChIP=
    unlist(lapply(1:nrow(region_in_factor),function(x){
      paste(unlist(strsplit(region_in_factor[x][[paste0(stat_in,'_motif')]],';')[[1]][which(unlist(lapply(strsplit(region_in_factor[x][[paste0(stat_in,'_motif')]],';')[[1]],function(mt){
        
        return(grepl(mt,region_in_factor[x]$motif_ChIP,ignore.case = T))
        
      })))]),collapse=';')
      
      
      
    }))
  
  }else(region_in_factor="No ChIP overlap")
  return(region_in_factor)
}
factor_olap<-function(tissue,factor_in,Ken_motif_folder,motif_locus_bed_dir=motif_locus_bed_dir,stage="embyro"){
  Ken_dNME=fread(paste0(Ken_motif_folder,tissue,'_OR_residual_dNME.csv'))
  Ken_dMML=fread(paste0(Ken_motif_folder,tissue,'_OR_residual_dMML.csv'))
  factor_in_dNME=factor_in[grepl(paste(gsub('.*_','',Ken_dNME$motif),collapse = "|"),factor_in$metadata,ignore.case=T)]
  factor_in_dMML=factor_in[grepl(paste(gsub('.*_','',Ken_dMML$motif),collapse = "|"),factor_in$metadata,ignore.case=T)]
  #Write the bedfile to output
  for(motif in gsub('.*_','',Ken_dNME$motif)){
    motif_ChIP=factor_in_dNME[grepl(motif,metadata,ignore.case = T)]
    if(nrow(motif_ChIP)>0){
      write.table(motif_ChIP,
                  paste0(motif_locus_bed_dir,'dNME/',tissue,'/mm10_',tissue,'_',stage,'_dNME_ChiPatlas_',gsub('::','_',motif),'.bed'),
                  col.names = F,row.names = F,quote=F)
    }
    
  }
  for(motif in gsub('.*_','',Ken_dMML$motif)){
    motif_ChIP=factor_in_dMML[grepl(motif,metadata,ignore.case = T)]
    if(nrow(motif_ChIP)>0){
      write.table(motif_ChIP,
                  paste0(motif_locus_bed_dir,'dMML/',tissue,'/mm10_',tissue,'_',stage,'_dMML_ChiPatlas_',gsub('::','_',motif),'.bed'),
                  col.names = F,row.names = F,quote=F)
    }
    
  }
  saveRDS(factor_in_dNME,paste0(ChiP_motif_dir,'factor_in_',stage,'_',tissue,'_all_dNME.rds'))
  saveRDS(factor_in_dMML,paste0(ChiP_motif_dir,'factor_in_',stage,'_',tissue,'_all_dMML.rds'))
  return(list(factor_in_dNME=factor_in_dNME,
              factor_in_dMML=factor_in_dMML))
  
}
ChIP_olap<-function(factor_in_dNME,factor_in_dMML,region_dNME,region_dMML){
  
  if(nrow(factor_in_dNME)>0){
    region_dNME_tissue_factor=ChIP_assignment(region_dNME,factor_in_dNME,stat_in="dNME")
  }else(region_dNME_tissue_factor="No significant factors")
  if(nrow(factor_in_dMML)>0){
    region_dMML_tissue_factor=ChIP_assignment(region_dMML,factor_in_dMML,stat_in="dMML")
  }else(region_dMML_tissue_factor="No significant factors")
  return(list(dNME_region=region_dNME_tissue_factor,
              dMML_region=region_dMML_tissue_factor))
  
  
}
