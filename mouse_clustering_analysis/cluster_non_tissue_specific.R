source('mainFunctions_sub.R')
cut <- 0.1
d=readRDS(UC_in_matrix_cluster_file)
# lapply(d,nrow)
# $EFP
# [1] 4876367
# $forebrain
# [1] 4876367
# $heart
# [1] 4876367
# $hindbrain
# [1] 4876367
# $limb
# [1] 4876367

# $liver
# [1] 4876367
# $midbrain
# [1] 4876367
aid <- sapply(names(d),function(i) {
    names(which(rowSums(d[[i]] > cut) > 0))
  })  
#lapply(aid,length)
# $EFP
# [1] 935740
# $forebrain
# [1] 1165702
# $heart
# [1] 994001
# $hindbrain
# [1] 1132357
# $limb
# [1] 970433
# $liver
# [1] 1093959
# $midbrain
# [1] 1113027
timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
for (seed in 1:10) {
  cat('Processing:',seed,'\n')
  cluster_d <- sapply(names(d),function(i) {
    #This is for one and only one
    #sid <- setdiff(aid[[i]],unlist(aid[names(aid)!=i]))
    sid=aid[[i]]
    i <- d[[i]]
    i <- i[sid,]
    i <- i[,colnames(i) %in% timeorder]
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
  dir_uc=paste0(dir_cluster_in_non_ts,'uc_',sub('\\.','',as.character(cut)),'/')
    ifelse(!dir.exists(file.path(dir_uc)), dir.create(file.path(dir_uc)), FALSE)
  saveRDS(cluster_d,file=paste0(dir_uc,'uc_',cut,'_',seed,'.rds'))
}




