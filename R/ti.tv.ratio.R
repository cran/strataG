#' @rdname ti.tv.ratio
#' @export ti.tv.ratio is.ti.tv
#' 
#' @title Transition / Transversion Ratio
#' @description Calculate transition/transversion ratio. Test substitution type of two bases.
#' 
#' @param dna.seq sequence of DNA.
#' @param b1,b2 two bases to be compared.
#' @param sub.type type of substitution to check: 'ti' for transition, 'tv' for transversion.
#' @return For \code{ti.tv.ratio}, a vector providing: \cr
#' \tabular{ll}{
#'   \code{Ti} \tab the number of transitions.\cr
#'   \code{Tv} \tab the number of transversions.\cr
#'   \code{Ti.Tv.ratio} \tab the transition/transversion ratio.\cr
#' }
#' For \code{is.ti.tv} a logical identifying whether the \code{b1} to \code{b2} is a \code{sub.type} substitution
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}

ti.tv.ratio <- function(dna.seq) {
  stopifnot.aligned(dna.seq)
  
  site.freqs <- base.freqs(dna.seq, c("a", "c", "g", "t"))$site.freqs
  ti.tv <- apply(site.freqs, 2, function(freqs) {
    freqs <- freqs[freqs > 0]
    if(length(freqs) < 2) {
      c(Ti = 0, Tv = 0)
    } else {
      pairs <- combn(names(freqs), 2)
      c(Ti = sum(apply(pairs, 2, function(bases) is.ti.tv(bases[1], bases[2], "ti"))),
        Tv = sum(apply(pairs, 2, function(bases) is.ti.tv(bases[1], bases[2], "tv")))
      )
    }
  })
  
  Ti = sum(ti.tv["Ti", ])
  Tv = sum(ti.tv["Tv", ])
  c(Ti = Ti, Tv = Tv, Ti.Tv.ratio = Ti / Tv)
}