af <- list.files('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/saver/c1')
for (f in af) {
  expr <- readRDS(paste0('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/data/saver/c1/',f))
  expr <- log2(expr + 1)
  expr <- expr[rowMeans(expr > 0.1) > 0.01,]
    expr <- expr[!grepl('^Rpl|^Rps',rownames(expr)),]
  m <- rowMeans(expr)
  lm <- log(m)
  var <- apply(expr,1,var)
  logvar <- log(var)
  hypervar_var <- resid(loess(var~lm))
  hypervar_logvar <- resid(loess(logvar~lm))
  res <- data.frame(mean=m,var=var,hypervar_var=hypervar_var,hypervar_logvar=hypervar_logvar)
  saveRDS(res,file=paste0('/home-4/zji4@jhu.edu/scratch/andy_ASE/hypervar/encmouse/saverres/c1/',f))
}

