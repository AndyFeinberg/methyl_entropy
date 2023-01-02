source('mainFunctions_sub.R')
UC_in=readRDS(UC_merge_file)
dnme=mclapply(UC_in,function(x) {x=x[,grepl('dNME',colnames(x))];colnames(x)=gsub('dNME-|-all','',colnames(x));return(x)},mc.cores=7)
dmml=mclapply(UC_in,function(x) {x=x[,grepl('dMML',colnames(x))];colnames(x)=gsub('dMML-|-all','',colnames(x));return(x)},mc.cores=7)
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
  #i <- i[rowSums(i) > 0,]
  i <- i[,colnames(i) %in% timeorder]
  i <- i[,order(match(colnames(i),timeorder))]
  i <- scalematrix(i)
  i <- i[complete.cases(i),]
})

dmmlcor <- dnmecor <- list()
for (n in names(d)) {
  dmmlcor[[n]] <- corfunc(dmml[[n]][rownames(d[[n]]),colnames(d[[n]])],d[[n]])
  dnmecor[[n]] <- corfunc(dnme[[n]][rownames(d[[n]]),colnames(d[[n]])],d[[n]])
}

saveRDS(dmmlcor,file=dmml_cor_file)
saveRDS(dnmecor,file=dnme_cor_file)


