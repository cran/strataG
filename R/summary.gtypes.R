#' @export summary.gtypes
#' 
#' @title Summarize gtypes Object
#' @description Generate a summary of a \code{gtypes} object.
#' 
#' @usage \method{summary}{gtypes}(object, ... )
#' 
#' @param object a \code{\link{gtypes}} object.
#' @param ... other arguments (ignored).
#' 
#' @return a list with the following elements:
#' \tabular{ll}{
#'   \code{all} \tab vector of number of samples, strata, and loci.\cr
#'   \code{locus.freq} \tab for haploid \code{gtypes}, a frequency table of haplotypes by strata, and for diploid.\cr
#'   \code{gtypes} \tab a list of strata containing a list of data.frames of allele frequencies for each loci.\cr
#'   \code{by.strta} \tab a by-strata data.frame summarizing haplotypes or loci.\cr
#'   \code{num.seq} \tab if haploid, the number of sequences, otherwise \code{NA}.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

summary.gtypes <- function(object, ...) {  
  stopifnot.gtypes(object)
  strata.names <- attr(object, "strata.names")
  
  locus.freq <- if(is.diploid(object)) {
    sapply(strata.split(object, remove.sequences = TRUE), allele.freqs, simplify = FALSE)
  } else hap.freqs(object)
  
  gtype.summary <- c(num.samples = nrow(object$genotypes), 
                     num.strata = length(unique(decode.strata(object))),
                     num.loci = length(attr(object, "locus.names"))
  )
  
  if(is.haploid(object)) {
    gtype.summary["seq.length"] <- ifelse(is.aligned.seq(object$sequences), length(object$sequences[[1]]), NA)
  }
  
  strata.summary <- t(sapply(strata.split(object, remove.sequences = TRUE), function(strata) {
    if(is.haploid(strata)) {
      c(num.samples = nrow(strata$genotypes), 
        num.haps = num.alleles(strata), 
        hap.div = haplotypic.diversity(strata),
        pct.unique.haps = pct.unique.haplotypes(strata)
      )
    } else {
      c(num.samples = nrow(strata$genotypes),
        mean.num.alleles = mean(num.alleles(strata)),
        mean.allelic.richness = mean(allelic.richness(strata)),
        mean.heterozygosity = mean(obsvd.het(strata))
      )
    }
  }))
  num.seq <- ifelse(is.haploid(object), length(object$sequences), NA)
  
  sum.g <- list(all = gtype.summary, locus.freq = locus.freq, by.strata = strata.summary, num.seq = num.seq)
  attr(sum.g, "description") <- attr(object, "description")
  class(sum.g) <- c("gtypeSummary", "list")
  sum.g
} 