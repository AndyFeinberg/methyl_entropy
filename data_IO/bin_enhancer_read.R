source('mainFunctions_sub.R')
# Read in mouse enhancer from Bin's paper--------------------------------------------------
bin_supp=fread('../downstream/input/mouse_analysis/enhancer_selection/s8C_enhancer_bin.csv',skip=1)
bin_supp_gr=GRanges(seqnames=bin_supp$chrom,ranges=IRanges(start=bin_supp$start,end=bin_supp$end))
mcols(bin_supp_gr)=bin_supp[,c(-1,-2,-3)]
colnames(mcols(bin_supp_gr))=gsub('\\...*','',colnames(mcols(bin_supp_gr)))
saveRDS(bin_supp_gr,bin_enhancer_rds)
enhancer_Bin=readRDS(bin_enhancer_rds)
#Generate bed files
export(sort(enhancer_Bin),bin_enhancer_bed)

