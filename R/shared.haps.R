#' @export shared.haps
#' 
#' @title Shared Haplotypes
#' @description Identify haplotypes shared between strata.
#' 
#' @param g a \code{\link{gtypes}} object.
#' 
#' @return a data.frame identifying the hapltoypes shared between each pair of strata.
#'   If no haplotypes are shared, returns \code{NA}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

shared.haps <- function(g) {
  df <- as.data.frame(g)
  strata <- sort(unique(df$strata))
  if(length(strata) < 2) stop("'g' has less than 2 strata")
  strata.pairs <- as.data.frame(t(combn(strata, 2)), stringsAsFactors = FALSE)
  colnames(strata.pairs) <- c("strata.1", "strata.2")
  
  strata.pairs$shared.haps <- sapply(1:nrow(strata.pairs), function(i) {
    pair.df <- subset(df, strata %in% unlist(strata.pairs[i, ]))
    freq.table <- table(df$haplotype, df$strata)
    shared <- apply(freq.table, 1, function(i) sum(i > 0) == 2)
    if(sum(shared) == 0) {
      NA
    } else {
      paste(rownames(freq.table)[shared], collapse = ", ")
    }
  })
  strata.pairs  
}