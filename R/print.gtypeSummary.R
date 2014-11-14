#' @export print.gtypeSummary
#' 
#' @title Print gtypeSummary
#' @description Print the output for a \code{gtypeSummary} object from \code{\link{summary.gtypes}}.
#' 
#' @usage \method{print}{gtypeSummary}(x, ... )
#' 
#' @param x list from summary.gtypes
#' @param ... other arguments (ignored)
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

print.gtypeSummary <- function(x, ... ) { 
  cat("\n") 
  desc <- attr(x, "description")
  if(!is.null(desc)) cat("--", desc, "--\n")
  cat("Samples:", x$all["num.samples"], "\n")
  cat("Strata:", x$all["num.strata"], "\n")
  loc.names <- if(is.table(x$locus.freq)) {
    "Haplotype"
  } else {
    names(x$locus.freq[[1]])
  }
  cat("Loci:", x$all["num.loci"], "-", loc.names, "\n")
  if(!is.na(x$num.seq)) cat("Number of sequences:", x$num.seq, "\n")
  cat("\n")
  print(x$by.strata)
  cat("\n")
  invisible(x)
}