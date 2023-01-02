source('mainFunctions_sub.R')
UC_in=readRDS(UC_merge_file)
dnme_all=mclapply(UC_in,function(x) {x=x[,grepl('dNME',colnames(x))];colnames(x)=gsub('dNME-|-all','',colnames(x));return(x)},mc.cores=7)
dmml_all=mclapply(UC_in,function(x) {x=x[,grepl('dMML',colnames(x))];colnames(x)=gsub('dMML-|-all','',colnames(x));return(x)},mc.cores=7)
d=mclapply(UC_in,function(x) {x=x[,grepl('UC',colnames(x))];colnames(x)=gsub('UC-|-all','',colnames(x));colnames(x)=sub('.*?-','',colnames(x));return(x)},mc.cores=7)

scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}

corfunc <- function(m1,m2,type='concordant') {
  if (type=='concordant') {
    rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
  } else {
    scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)            
  }
}

timeorder <- sapply(1:20,function(i) paste0('E',i,'.5-E',i+1,'.5'))
d <- sapply(d,function(i) {
  #i <- i[rowSums(i > 0.1) > 0,]
  i <- i[,colnames(i) %in% timeorder]
  i <- i[,order(match(colnames(i),timeorder))]
  i <- scalematrix(i)
  i <- i[complete.cases(i),]
})

dmmlcor <- dnmecor <- list()
 library(combinat)
for (n in names(d)) {
  dmml=dmml_all[[n]][rownames(d[[n]]),colnames(d[[n]])]
  dnme=dnme_all[[n]][rownames(d[[n]]),colnames(d[[n]])]
  mat <- do.call(rbind,combinat::permn(1:ncol(d[[n]])))
  iden <- apply(mat,1,function(i) {
    mean(i==sort(i))
  })
  mat <- mat[iden!=1,]
  set.seed(12345)
  mat <- mat[sample(1:nrow(mat),20),]
  sampid <- t(mat)
  print(sum(apply(sampid,2,function(i) {
    mean(i==sort(i))
  })==1))
  dmmlcor[[n]] <- sapply(1:20,function(i) corfunc(dmml,d[[n]][,sampid[,i]]))
  dnmecor[[n]] <- sapply(1:20,function(i) corfunc(dnme,d[[n]][,sampid[,i]]))
}
dir_in_cor_Jason='../downstream/input/mouse_analysis/correlation_analysis/all_regions/'
saveRDS(dmmlcor,file=paste0(dir_in_cor_Jason,'fullpermudmmlcor_YQ.rds'))
saveRDS(dnmecor,file=paste0(dir_in_cor_Jason,'fullpermudnmecor_YQ.rds'))



