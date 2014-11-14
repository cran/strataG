#' @rdname pop.diff.test
#' @importFrom parallel mclapply 
#' @importFrom swfscMisc one.arg

overall.test <- function(g, stats = "all", nrep = 1000, keep.null = FALSE, num.cores = 1, quietly = FALSE, write.output = FALSE, ...) {  
  # check gtypes
  g <- check.gtypes(g)
  
  # check stats
  stats <- stat.list(g, stats)
  
  # check replicates
  if(!is.numeric(nrep) & length(nrep) != 1) stop("'nrep' must be a single-element numeric vector")
  if(nrep < 1) keep.null = FALSE
  
  if(!quietly) cat(format(Sys.time()), ":", attr(g, "description"), 
                   ": Running", nrep, "permutations for", length(stats), 
                   ifelse(length(stats) == 1, "analysis...", "stats..."), "\n"
  )
  
  # calculate list of observed values for each population structure function
  obs <- lapply(stats, function(f) if(one.arg(f)) f(g) else f(g, ...))
  names(obs) <- sapply(obs, function(x) x$stat.name)
  
  # calculate matrix of null distributions (rows are reps, columns are metrics)
  opt <- options(check.gtypes = FALSE)
  null.dist <- if(nrep > 0) {
    null.g <- g
    # **** change this to create matrix of random strata assignments and mclapply through that sending
    # **** random strata vector and locus data to each C function directly
    do.call(rbind, mclapply(1:nrep, function(i) {
      null.g$genotypes[, "strata"] <- sample(g$genotypes[, "strata"])
      sapply(stats, function(f) if(one.arg(f)) f(null.g)$estimate else f(null.g, ...)$estimate)
    }, mc.cores = num.cores))
  } else rbind(sapply(stats, function(f) NA))  

  colnames(null.dist) <- names(obs)
  options(opt)
  
  # calculate vector of p-values
  pval <-sapply(names(obs), function(stat) p.val(obs[[stat]]$estimate, null.dist[, stat]))
  names(pval) <- names(obs)
  
  null.dist <- if(keep.null) null.dist else NULL
  
  # format matrix of estimates and p-values
  result <- t(sapply(names(obs), function(stat) c(obs[[stat]]$estimate, pval[stat])))
  colnames(result) <- c("estimate", "p.val")
  
  # collect strata frequencies to named vector (assuming that frequencies from first analysis are same for all stats)
  strata.freq <- as.numeric(obs[[1]]$strata.freq)
  names(strata.freq) <- names(obs[[1]]$strata.freq)
  
  if(!quietly) {
    sf <- cbind(strata.freq)
    colnames(sf) <- "N"
    cat("\n")
    print(sf)
    cat("\nPopulation structure results:\n")
    print(result)
    cat("\n")
  }
  
  if(write.output) {
    out.file <- gsub("[[:punct:]]", ".", paste(attr(g, "description"), "permutation test results.csv"))
    write.csv(result, out.file)
  }
  
  options(opt)
  invisible(list(strata.freq = strata.freq, result = result, null.dist = null.dist))
}
