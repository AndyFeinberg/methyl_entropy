motifbreakR_parallel<-function (snpList, pwmList, threshold = 0.85, filterp = FALSE, 
          method = "default", show.neutral = FALSE, verbose = FALSE, 
          bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25), legacy.score = TRUE, 
          BPPARAM = bpparam()) 
{
  if (.Platform$OS.type == "windows" && inherits(BPPARAM, 
                                                 "MulticoreParam")) {
    warning(paste0("Serial evaluation under effect, to achive parallel evaluation under\n", 
                   "Windows, please supply an alternative BPPARAM"))
  }
  cores <- bpnworkers(BPPARAM)
  num.snps <- length(snpList)
  if (num.snps < cores) {
    cores <- num.snps
  }
  if (is(BPPARAM, "SnowParam")) {
    bpstart(BPPARAM)
    cl <- bpbackend(BPPARAM)
    clusterEvalQ(cl, library("MotifDb"))
  }
  genome.package <- attributes(snpList)$genome.package
  if (requireNamespace(eval(genome.package), quietly = TRUE, 
                       character.only = TRUE)) {
    genome.bsgenome <- eval(parse(text = paste(genome.package, 
                                               genome.package, sep = "::")))
  }
  else {
    stop(paste0(eval(genome.package), " is the genome selected for this snp list and \n", 
                "is not present on your environment. Please load it and try again."))
  }
  pwms <- motifbreakR:::preparePWM(pwmList = pwmList, filterp = filterp, 
                     scoreThresh = threshold, bkg = bkg, method = method)
  k <- max(sapply(pwms$pwmList, ncol))
  snpList <- motifbreakR:::prepareVariants(fsnplist = snpList, genome.bsgenome = genome.bsgenome, 
                             max.pwm.width = k, legacy = legacy.score)
  snpList_cores <- split(as.list(rep(names(snpList), times = cores)), 
                         1:cores)
  for (splitr in seq_along(snpList)) {
    splitcores <- sapply(suppressWarnings(split(snpList[[splitr]], 
                                                1:cores)), list)
    for (splitcore in seq_along(snpList_cores)) {
      snpList_cores[[splitcore]][[splitr]] <- splitcores[[splitcore]]
      names(snpList_cores[[splitcore]])[splitr] <- names(snpList)[splitr]
    }
  }
  snpList <- snpList_cores
  rm(snpList_cores)
  cat('score SNPlist\n')
  #Changed
  #Check if modification happened here
  x <- bplapply(snpList, motifbreakR:::scoreSnpList, pwmList = pwms$pwmList, 
              threshold = pwms$pwmThreshold, pwmList.pc = pwms$pwmListPseudoCount, 
              pwmRanges = pwms$pwmRange, method = method, bkg = bkg, 
              show.neutral = show.neutral, verbose = ifelse(cores == 
                                                              1, verbose, FALSE), genome.bsgenome = genome.bsgenome, 
              filterp = filterp,BPPARAM=BPPARAM)
  if (inherits(x, "try-error")) {
    if (is(BPPARAM, "SnowParam")) {
      bpstop(BPPARAM)
    }
    stop(attributes(x)$condition)
  }
  if (is(BPPARAM, "SnowParam")) {
    bpstop(BPPARAM)
  }
  drops <- sapply(x, is.null)
  x <- x[!drops]
  pwmList <- pwms$pwmList
  pwmList@listData <- lapply(pwms$pwmList, function(pwm) {
    pwm <- pwm[c("A", "C", "G", "T"), 
               ]
    return(pwm)
  })
  pwmList.pc <- lapply(pwms$pwmListPseudoCount, function(pwm) {
    pwm <- pwm[c("A", "C", "G", "T"), 
               ]
    return(pwm)
  })
  if (length(x) > 1) {
    x <- unlist(GRangesList(unname(x)))
    snpList <- unlist(GRangesList(lapply(snpList, `[[`, 
                                         "fsnplist")), use.names = F)
    x <- x[order(match(names(x), names(snpList)), x$geneSymbol), 
           ]
    attributes(x)$genome.package <- genome.package
    attributes(x)$motifs <- pwmList[mcols(pwmList)$providerId %in% 
                                      unique(x$providerId) & mcols(pwmList)$providerName %in% 
                                      unique(x$providerName), ]
    attributes(x)$scoremotifs <- pwmList.pc[names(attributes(x)$motifs)]
  }
  else {
    if (length(x) == 1L) {
      x <- x[[1]]
      attributes(x)$genome.package <- genome.package
      attributes(x)$motifs <- pwmList[mcols(pwmList)$providerId %in% 
                                        unique(x$providerId) & mcols(pwmList)$providerName %in% 
                                        unique(x$providerName), ]
      attributes(x)$scoremotifs <- pwmList.pc[names(attributes(x)$motifs)]
    }
    else {
      warning("No SNP/Motif Interactions reached threshold")
      x <- NULL
    }
  }
  if (verbose && cores > 1) {
    if (is.null(x)) {
      message(paste("reached end of SNPs list length =", 
                    num.snps, "with 0 potentially disruptive matches to", 
                    length(unique(x$geneSymbol)), "of", length(pwmList), 
                    "motifs."))
    }
    else {
      message(paste("reached end of SNPs list length =", 
                    num.snps, "with", length(x), "potentially disruptive matches to", 
                    length(unique(x$geneSymbol)), "of", length(pwmList), 
                    "motifs."))
    }
  }
  return(x)
}