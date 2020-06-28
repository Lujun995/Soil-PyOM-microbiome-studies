GMPR<-function (comm, intersect.no = 4, ct.min = 2, verbose = FALSE) {
  # Computes the GMPR size factor
  #
  # Args:
  #   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
  #   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
  #   ct.min: the minimum number of counts required to calculate ratios ct.min = 5 has better results
  
  #
  # Returns:
  #   a list that contains:
  #      gmpr: GMPR size factors for all samples; Samples with distinct sets of features will be output as NA.
  #      nss:   number of samples with significant sharing (> intersect.no) including itself
  
  # mask counts < ct.min
  comm[comm < ct.min] <- 0
  
  if (is.null(colnames(comm))) {
    colnames(comm) <- paste0('S', 1:ncol(comm))
  }
  
  if (verbose == TRUE)
    cat('Begin GMPR size factor calculation ...\n')
  
  comm.no <- numeric(ncol(comm))
  gmpr <- sapply(1:ncol(comm),  function(i) {
    
    if (i %% 50 == 0) {
      if (verbose == TRUE)
        cat(i, '\n')
    }
    x <- comm[, i]
    # Compute the pairwise ratio
    pr <- x / comm
    # Handling of the NA, NaN, Inf
    pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
    # Counting the number of non-NA, NaN, Inf
    incl.no <- colSums(!is.na(pr))
    # Calculate the median of PR
    pr.median <- colMedians(pr, na.rm=TRUE)
    # Record the number of samples used for calculating the GMPR
    comm.no[i] <<- sum(incl.no >= intersect.no)
    # Geometric mean of PR median
    if (comm.no[i] > 1) {
      return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
    } else {
      return(NA)
    }
  }
  )
  
  if (sum(is.na(gmpr))) {
    warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'), 
                   '\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
                   'For these samples, their size factors are set to be NA! \n', 
                   'You may consider removing these samples since they are potentially outliers or negative controls!\n',
                   'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
  }
  if (verbose == TRUE) {
    cat('Completed!\n')
    cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
  }
  
  attr(gmpr, 'NSS') <- comm.no
  names(gmpr) <- colnames(comm)
  return(gmpr * median(colSums(comm)))
}


perm_fdr_adj<-function (F0, Fp) {
  ord <- order(F0, decreasing = T)
  F0 <- F0[ord]
  perm.no <- ncol(Fp)
  Fp <- as.vector(Fp)
  Fp <- Fp[!is.na(Fp)]
  Fp <- sort(c(Fp, F0), decreasing = F)
  n <- length(Fp)
  m <- length(F0)
  FPN <- (n + 1) - match(F0, Fp) - 1:m
  p.adj.fdr <- FPN / perm.no / (1:m)
  #p.adj.fdr <- sapply(F0, function(x) sum(Fp >= 
  #x, na.rm=TRUE) / perm.no)/(1:length(F0))
  p.adj.fdr <- pmin(1, rev(cummin(rev(p.adj.fdr))))[order(ord)]
}

perm_fwer_adj<-function (F0, Fp) {
  ord <- order(F0, decreasing = T)
  F0 <- F0[ord]
  col.max <- colMaxs(Fp, na.rm=TRUE)
  p.adj.fwer <- sapply(F0, function(x) mean(col.max >= x))[order(ord)]
}

na.pad<-function (vec, ind) {
  vec0 <- numeric(length(ind))
  vec0[!ind] <- vec
  vec0[ind] <- NA
  vec0
}

permute_differential_analysis<-function (meta.dat, comm, grp.name, adj.name = NULL, size.factor = NULL, 
            transform = 'arcsqrt', weights = NULL, strata = NULL,  perm.no = 999, 
            stage.no = 1, stage.pv = 0.05, stage.max.pct = 0.20, verbose = TRUE) {
  
  # Args:
  #   meta.dat: a data frame containing the sample information
  #   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
  #   size.factor:  a numeric vector of the library sizes; if NULL, GMPR size factors will be used
  #   weights: a vector of the weights; if null, the data will be weighted by size factor
  #   grp.name: a character, variable of interest; it could be numeric or categorical; Should be in meta.dat
  #   adj.name: a character vector, variable(s) to be adjusted; they could be numeric or categorical; Should be in meta.dat
  #   strata:  a factor indicating the permutation strata; permutation will be confined to each stratum 
  #   perm.no: the number of permutations; If the FDR/FWER-adjusted p values are the major interest, 
  #           perm.no could be set to 50 to reduce computation
  #   stage.no: the number of stages if multiple-stage normalization stategy is used
  #   stage.pv: the raw p value cutoff below which the features will be excluded for calculating the size factor
  #   stage.max.pct: the maximum percentage of features that will be excluded
  #   verbose: whether the trace information should be printed out
  
  #
  # Returns:
  #   a list that contains:
  #      call: the call
  #      R2: a vector of percent explained variance for 
  #      p.value: the raw p-values based on permutations
  #      p.adj.fdr: permutation-based FDR-adjusted p.value
  #      p.adj.fwer: permutation-based FWER-adjusted p.value
  #      size.factor: the size.factor used
  #      weights: the weights used
  
  this.call = match.call()
  
  if (is.null(size.factor)) {
    size.factor <- GMPR(comm)
  } else {
  }
  
  n <- ncol(comm)
  row.names <- rownames(comm)
  
  for (i in 1:stage.no) {
    if (verbose == TRUE)
      cat('Stage ', i, '...\n')
    if (is.null(weights)) {
      W <- size.factor
    } else {
      W <- weights*size.factor
    }
    
    W <- sqrt(W)
    
    Y <- t(t(comm) / size.factor)
    
    if (transform == 'arcsqrt') {
      Y[Y <= 0] <- 0
      Y[Y >= 1] <- 1
      Y <- asin(sqrt(Y))
    } 
    
    Y <- W * Y
    # Covariate space (including intercept)
    
    if (is.null(adj.name)) {
      M0 <- model.matrix(~ 1, meta.dat) 
    } else {
      df0 <- meta.dat[, c(adj.name), drop = FALSE]
      M0 <- model.matrix( ~ ., df0) 
    }
    
    M0 <- W * M0
    
    # Remove covariate effects
    Y <- t(resid(lm(as.formula(paste('t(Y) ~ M0 - 1')), meta.dat)))
    
    if (!is.null(strata)) {
      strata <- factor(strata)
    }
    
    # Residual space after adjusting covariate
    df1 <- meta.dat[, c(grp.name), drop = FALSE]
    M1 <- model.matrix( ~ . - 1, df1) 
    M1 <- W * M1
    M1 <- as.matrix(resid(lm(M1 ~ M0 - 1)))
    
    # QR decompostion
    qrX0 <- qr(M0, tol = 1e-07)
    Q0 <- qr.Q(qrX0)
    Q0 <- Q0[, 1:qrX0$rank, drop = FALSE]
    
    qrX1 <- qr(M1, tol = 1e-07)
    Q1 <- qr.Q(qrX1)
    Q1 <- Q1[, 1:qrX1$rank, drop = FALSE] 
    
    TSS <- rowSums(Y^2)
    MSS1 <- rowSums((Y %*% Q1)^2)
    
    # Scaled F-stat
    F0 <- MSS1 /  (TSS - MSS1)  
    R2 <- MSS1 / TSS
    
    perm.ind <- vegan:::getPermuteMatrix(perm.no, n, strata = strata)
    perm.no <- nrow(perm.ind)
    
    Fp <- sapply(1:perm.no, function(i) {
      if (verbose) {
        if (i %% 100 == 0) cat('.')
      }
      Q1p <- Q1[perm.ind[i, ], , drop = FALSE]
      MSS1p <- rowSums((Y %*% Q1p)^2)
      MSS1p /  (TSS - MSS1p) 
    })
    if (verbose) {
      cat('\n')
    }
    
    if (mean(is.na(F0)) >= 0.1) {
      warning('More than 10% observed F stats have NA! Please check! \n')
    }
    
    if (mean(is.na(Fp)) >= 0.1) {
      warning('More than 10% permuted F stats have NA! Please check! \n')
    }
    
    na.ind <- is.na(F0)
    F0 <- F0[!na.ind]
    Fp <- Fp[!na.ind, ]
    
    p.raw <- cbind(Fp >= F0, 1)
    p.raw <- rowMeans(p.raw)
    
    if (i == stage.no) {
      break
    } else {
      # recalculating the size factor
      if (mean(p.raw <= stage.pv) > stage.max.pct) {
        ind <- p.raw > quantile(p.raw, stage.max.pct)
      } else {
        ind <- p.raw > stage.pv
      }
      size.factor <- GMPR(comm[ind, ])
    }
    
  }
  
  # scaled F stat
  p.adj.fdr <- perm_fdr_adj(F0, Fp)
  p.adj.fwer <- perm_fwer_adj(F0, Fp)
  
  #Fp <-  1 - (apply(Fp, 1, rank) - 1) / ncol(Fp)
  #Fp <- t(Fp)
  #p.adj.fdr <- perm_fdr_adj(-p.raw, -Fp)
  #p.adj.fwer <- perm_fwer_adj(-p.raw, -Fp)
  
  p.raw <- na.pad(p.raw, na.ind)
  p.adj.fdr <- na.pad(p.adj.fdr, na.ind)
  p.adj.fwer <- na.pad(p.adj.fwer, na.ind)
  
  names(p.raw) <- names(p.adj.fdr) <- names(p.adj.fwer) <- row.names
  
  if (verbose) cat('Completed!\n')
  return(list(call = this.call, R2 = R2,  p.raw = p.raw, p.adj.fdr = p.adj.fdr, p.adj.fwer = p.adj.fwer,
              size.factor = size.factor, weights = weights))
  
}


