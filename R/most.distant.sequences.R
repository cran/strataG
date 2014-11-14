#' @export most.distant.sequences
#' 
#' @title Most Distant Sequences
#' @description Finds the set of sequences that have the greatest mean pairwise distance and 
#'   smallest variance of pairwise distances.
#' 
#' @param x a \code{\link{gtypes}} object with aligned sequences, or a list of aligned DNA sequences.
#' @param num.seqs number of sequences to return.
#' @param model a character string specifying the evolutionary model to be used. 
#'   See \link{dist.dna} for more information.
#' @param pairwise.deletion a logical indicating whether to delete sites with missing data. 
#'   See \link{dist.dna} for more information.
#' 
#' @return a vector of the sequence names that are have the greatest mean pairwise distance and 
#'   smallest variance of pairwise distances. The names are returned in order from most to least distant.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

most.distant.sequences <- function(x, num.seqs = NULL, model = "raw", pairwise.deletion = TRUE) {  
  if(is.gtypes(x, FALSE)) x <- decode.sequences(x)
  stopifnot.aligned(x)
  if(is.null(num.seqs)) num.seqs <- length(x)
  
  # calculate distance between sequences
  seq.dist <- as.matrix(dist.dna(as.DNAbin(x), model = model, pairwise.deletion = pairwise.deletion))
  
  # convert distances to coordinates
  opt <- options(warn = -1)
  seq.cmd <- cmdscale(seq.dist, k = length(x) - 1)
  options(opt)
  
  # normalize coordinates to have mean of 0
  seq.cmd <- t(t(seq.cmd) - apply(seq.cmd, 2, function(p) mean(range(p, na.rm = TRUE))))
  
  # calculate euclidean distance to mean
  euc.dist <- sort(apply(seq.cmd, 1, function(p) sqrt(sum(p ^ 2, na.rm = TRUE))), decreasing = TRUE)
  
  # get the most distant 5 sequences
  ids <- names(euc.dist)[1]
  for(i in 1:(num.seqs - 1)) {      
    # add sequence with greatest mean distance and smallest variance
    other.seqs <- setdiff(rownames(seq.dist), ids)
    mean.var <- data.frame(t(sapply(other.seqs, function(hap) {
      m <- mean(seq.dist[hap, ids])
      v <- var(seq.dist[hap, ids])
      c(mean = m, var = v)
    })))
    mean.var <- cbind(mean.var, euc.dist = euc.dist[rownames(mean.var)])
    mean.var <- mean.var[order(mean.var$mean, -mean.var$var, mean.var$euc.dist, decreasing = TRUE), ]
    ids <- c(ids, rownames(mean.var)[1])
  }
  ids
}