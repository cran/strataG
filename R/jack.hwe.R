#' @rdname jack.hwe
#' @export jack.hwe jack.influential plot.jack.influential
#' 
#' @title Hardy-Weinberg Equlibrium Jackknife
#' @description Test influence of samples on Hardy-Weinberg equilibrium via jackknife.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param exclude.num Number of samples to exclude at a time.
#' @param min.hwe.samples minimum samples needed to calculate HWE.
#' @param show.progress logical. Show progress of jackknife?
#' @param ... other arguments to be passed to GENEPOP.
#' @param jack.result result from run of \code{jack.hwe}.
#' @param alpha critical value to determine if exclusion is "influential".
#' @param x result from a call to \code{jack.influential}.
#' @param main main title for influential sample plots from \code{plot.jack.influential}.
#' 
#' @details \tabular{ll}{
#'   \code{jack.hwe} \tab performs a HWE jackknife where all combinations of \code{exclude.num} samples 
#'     are left out and HWE is recalculated.\cr
#'   \code{jack.influential} \tab calculates odds.ratios between jackknife HWE and observed HWE and identifies
#'     "influential" samples. Samples are "influential" if the observed HWE p-value is < \code{alpha},
#'     but is > \code{alpha} when the samples are not present.\cr
#'   \code{plot.jack.influential} \tab creates a cumulative frequency plot of all odds-ratios from \code{jack.influential}. 
#'     A vertical dashed line marks the smallest influential exclusion.\cr
#' }
#' 
#' @return \code{jack.hwe} returns a list with:
#' \item{obs}{a named vector of HWE p-values for each locus.}
#' \item{jack}{a \code{data.frame} of HWE p-values where each row is an exclusion and columns are loci.}
#' \item{gtypes}{the original \code{gtypes} object.}\cr
#' \code{jack.influential} returns a list with:
#' \item{influential}{a \code{data.frame} of influential exclusions.}
#' \item{allele.freqs}{a \code{data.frame} listing the allele frequencies of influential exclusions.}
#' \item{odds.ratio}{a \code{matrix} of odds ratios between exclusions (rows) and loci (columns).}
#' 
#' @references Morin, P.A., R.G. LeDuc, F.I. Archer, K.K. Martien, R. Huebinger, J.W. Bickham,
#'   and B.L. Taylor. 2009. Significant deviations from Hardy-Weinberg equilibirum caused
#'   by low levesl of microsatellite genotyping errors. Molecular Ecology Resources 9:498-504.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

jack.hwe <- function(g, exclude.num = 1, min.hwe.samples = 5, show.progress = TRUE, ...) {  
  stopifnot.gtypes(g, "diploid")
  
  if((nrow(g$genotypes) - exclude.num) < min.hwe.samples){
    stop(paste("'exclude.num' or 'min.HWE.samples' is too large to analyze this data"))
  }
  
  warn <- getOption("warn")
  opt <- options(warn = -1, check.gtypes = FALSE)
  
  # Setup array of samples to exclude
  exclude.arr <- combn(rownames(g$genotypes), exclude.num)
  
  # Calculate observed HWE and run jackknife 
  start.time <- Sys.time()
  nsteps <- ncol(exclude.arr) + 1 
  if(show.progress) cat(format(start.time, "%Y-%m-%d %H:%M:%S"), "|", 1, "/", nsteps, "\n")
  obs <- hwe.genepop(g, show.output = FALSE, ...)
  jack <- sapply(1:ncol(exclude.arr), function(i) {           
    if(show.progress) {
      now <- Sys.time()    
      avg.time <- difftime(now, start.time, units = "secs") / (i + 1)
      est.complete.time <- now + (avg.time * (nsteps - i))
      cat(format(now, "%Y-%m-%d %H:%M:%S"), "|", i + 1, "/", nsteps, "| ETA =", 
          format(est.complete.time, "%Y-%m-%d %H:%M:%S\n")
      )
    }
    to.keep <- setdiff(rownames(g$genotypes), exclude.arr[, i])
    jack.gtypes <- subset(g, ids = to.keep)
    hwe.genepop(jack.gtypes, show.output = FALSE, ...)
  })

  exclude.vec <- apply(exclude.arr, 2, function(x) paste(x, collapse = ", "))
  options(warn = warn, check.gtypes = T)
  result <- list(obs = obs, jack = data.frame(excluded = exclude.vec, t(jack), stringsAsFactors = FALSE), gtypes = g)
  class(result) <- c("jack.hwe.result", class(result))
  options(opt)
  result
}