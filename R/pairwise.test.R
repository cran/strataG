#' @rdname pop.diff.test

pairwise.test <- function(g, stats = "all", nrep = 1000, keep.null = FALSE, num.cores = 1, quietly = FALSE, write.output = FALSE, ...) {
  g <- check.gtypes(g)
  
  # check stats
  stats <- stat.list(g, stats)
  
  # check replicates
  if(!is.numeric(nrep) & length(nrep) != 1) stop("'nrep' must be a single-element numeric vector")
  if(nrep < 1) keep.null = FALSE
  
  # create strata pairs
  strata <- attr(g, "strata.names")[-1]
  if(length(strata) < 2) stop("'g' has less than 2 strata")
  strata.pairs <- as.data.frame(t(combn(strata, 2)), stringsAsFactors = FALSE)
  colnames(strata.pairs) <- c("strata.1", "strata.2")
  
  if(!quietly) cat("\nRunning", nrep, "permutations for", length(stats), 
                   ifelse(length(stats) == 1, "analysis...", "stats..."), "\n"
  )
  
  opt <- options(check.gtypes = FALSE)
  
  # run permutation test on all pairwise gtypes subsets
  pair.list <- lapply(1:nrow(strata.pairs), function(i) {
    pair <- unlist(strata.pairs[i, ])
    if(!quietly) {
      pair.txt <- paste(attr(g, "description"), ":", paste(pair, collapse = " v. "))
      cat(format(Sys.time()), ":", pair.txt, "\n")
    }
    overall.test(g = subset(g, strata = pair), stats = stats, 
      nrep = nrep, keep.null = keep.null, quietly = TRUE, write.output = FALSE, 
      num.cores = num.cores, ...
    )
  })
  
  # compile results in 'pair.list' into a data.frame
  result <- do.call(rbind, lapply(pair.list, function(pair) {
    result.vec <- as.vector(t(pair$result))
    names(result.vec) <- paste(rep(rownames(pair$result), each = 2), c("", ".p.val"), sep = "") 
    result.vec <- rbind(result.vec)
    s1 <- names(pair$strata.freq)[1]
    s2 <- names(pair$strata.freq)[2]
    n1 <- pair$strata.freq[1]
    n2 <- pair$strata.freq[2]
    strata.1 <- paste(s1, " (", n1, ")", sep = "")
    strata.2 <- paste(s2, " (", n2, ")", sep = "")
    pair.label <- data.frame(pair.label = paste(strata.1, " v. ", strata.2, sep = ""))
    cbind(data.frame(pair.label = pair.label, strata.1 = s1, strata.2 = s2, n.1 = n1, n.2 = n2), result.vec)   
  }))
  rownames(result) <- NULL
 
  if(!quietly) {
    cat("\nPopulation structure results:\n")
    print(result[, c(1, 6:ncol(result))])
    cat("\n")
  }
  
  # create pairwise matrices - upper right is estimate, lower left is p-value
  stat.cols <- seq(6, ncol(result), 2)
  pair.mat <- lapply(stat.cols, function(i) {
    mat <- matrix(nrow = length(strata), ncol = length(strata), dimnames = list(strata, strata))
    for(j in 1:nrow(result)) {
      strata.1 <- as.character(result$strata.1[j])
      strata.2 <- as.character(result$strata.2[j])
      mat[strata.1, strata.2] <- result[j, i]
      mat[strata.2, strata.1] <- result[j, i + 1]
    }
    mat
  })
  names(pair.mat) <- colnames(result)[stat.cols]
  
  # compile null distributions into list of matrices
  null.dist <- if(keep.null) {
    null.mat <- lapply(pair.list, function(pair) pair$null.dist)
    names(null.mat) <- result$pair.label
    null.mat
  } else NULL
  
  if(write.output) {
    for(stat in names(pair.mat)) {
      out.file <- gsub("[[:punct:]]", ".", paste(attr(g, "description"), stat, "pairwise matrix.csv"))
      write.csv(pair.mat[[stat]], out.file)
    }
  }
  
  options(opt)
  invisible(list(result = result, pair.mat = pair.mat, null.dist = null.dist))
}