#' @export sequence.summary
#' 
#' @title Sequence Summary
#' @description Summaries for each sequence.
#' 
#' @param x a haploid \code{\link{gtypes}} object with sequences or a list of DNA sequences.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

sequence.summary <- function(x) {
  if(is.gtypes(x, FALSE)) x <- decode.sequences(x)
  if(!is.dna.seq(x, FALSE)) stop("'x' must contain valid sequences")
  
  data.frame(t(sapply(x, function(this.seq) {
    seq.rle <- rle(this.seq)
    start <- ifelse(seq.rle$values[1] == "-", seq.rle$lengths[1] + 1, 1)
    n <- length(seq.rle$values)
    sum.n <- sum(seq.rle$lengths)
    end <- ifelse(seq.rle$values[n] == "-", sum.n - seq.rle$lengths[n], sum.n)
    num.ns <- sum(this.seq == "n")
    num.indels <- sum(this.seq[start:end] == "-")
    c(start = start, end = end, length = end - start + 1, num.ns = num.ns, num.indels = num.indels)
  })))
}
