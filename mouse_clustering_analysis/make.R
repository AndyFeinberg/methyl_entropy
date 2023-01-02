library(RColorBrewer)
library(pheatmap)
setwd('../')
source('mainFunctions_sub.R')
cut <- as.numeric(commandArgs(trailingOnly = T))
for (seed in 1:10) {
  #see fulluc.R
  d <- readRDS(UC_in_matrix_cluster_file)
  timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
  
  aid <- sapply(names(d),function(i) {
    names(which(rowSums(d[[i]] > cut) > 0))
  })  
  
  d <- sapply(names(d),function(i) {
    #This is for one and only one
    sid <- setdiff(aid[[i]],unlist(aid[names(aid)!=i]))
    i <- d[[i]]
    i <- i[sid,]
    i <- i[,colnames(i) %in% timeorder]
    scalematrix <- function(data) {
      cm <- rowMeans(data)
      csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
      (data - cm) / csd
    }
    i <- scalematrix(i)
    i <- i[complete.cases(i),]
    set.seed(seed)
    clu <- kmeans(i,10,iter.max = 10000)$cluster
    n <- names(clu)
    clum <- rowsum(i,clu)/as.vector(table(clu))
    maxp <- apply(clum,1,function(i) {
      i[i < 0] <- 0
      i <- i-min(i)
      i <- i/max(i)
      sum(i*c(1:length(i)))/sum(i)
    })
    clu <- rank(maxp)[clu]
    names(clu) <- n
    clu
  })
  
  saveRDS(d,file=paste0('../downstream/input/mouse_analysis/clustering/tissue_specific/uc_',as.character(gsub('\\.','',cut)),'/uc_',cut,'_',seed,'.rds'))
}




