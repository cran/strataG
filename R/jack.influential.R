#' @rdname jack.hwe
#' @importFrom swfscMisc odds

jack.influential <- function(jack.result, alpha = 0.05) {  
  if(!"jack.hwe.result" %in% class(jack.result)) stop("'jack.result' is not from 'jack.hwe'")
  warn <- getOption("warn")
  opt <- options(warn = -1)
  
  obs <- jack.result$obs
  exclude.vec <- jack.result$jack$excluded
  jack <- t(as.matrix(jack.result$jack[, -1]))
  
  # Find influential samples and calculate odds and odds ratios
  obs.arr <- matrix(obs, length(obs), length(exclude.vec))
  which.infl <- which((obs.arr <= alpha) & (jack > alpha), arr.ind = TRUE)
  influential <- if(nrow(which.infl) > 0) {
    infl <- data.frame(excluded = exclude.vec[which.infl[, "col"]], stringsAsFactors = FALSE) 
    infl$locus <- names(obs)[which.infl[, "row"]]
    infl$obs.pval <- obs[which.infl[, "row"]]
    infl$jack.pval <- apply(which.infl, 1, function(i) jack[i[1], i[2]])
    infl$obs.odds <- odds(infl$obs.pval)
    infl$jack.odds <- odds(infl$jack.pval)
    infl$odds.ratio <- infl$jack.odds / infl$obs.odds
    infl <- infl[order(infl$odds.ratio, decreasing = TRUE), ] 
    rownames(infl) <- 1:nrow(infl) 
    infl  
  } else NULL
  
  # Calculate allele frequencies of influential samples 
  allele.freqs <- if(!is.null(influential)) {
    samples.loci <- do.call(rbind, lapply(1:nrow(influential), function(i) {
      samples <- unlist(strsplit(influential$excluded[i], ", "))
      data.frame(id = samples, locus = influential$locus[i], stringsAsFactors = FALSE)
    }))
    formatted.freqs <- allele.freq.format(unique(samples.loci), jack.result$gtypes)
    rownames(formatted.freqs) <- 1:nrow(formatted.freqs)
    colnames(formatted.freqs)[1] <- "id"
    formatted.freqs
  } else NULL
  
  # Calculate full odds ratio matrix
  odds.ratio <- t(odds(jack) / odds(obs.arr))
  rownames(odds.ratio) <- exclude.vec
  colnames(odds.ratio) <- names(obs)
  
  result <- list(influential = influential, allele.freqs = allele.freqs, odds.ratio = odds.ratio)
  class(result) <- c("jack.influential", class(result))
  options(opt)
  result
}