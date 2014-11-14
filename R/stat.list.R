#' @rdname pop.diff.test

stat.list <- function(g, stats = "all") {
  if(is.character(stats)) {
    # Check that stats specified are valid and appropriate for ploidy
    stats <- tolower(stats)
    haploid.stats <- c("fst", "phist", "chi2")
    diploid.stats <- c("fst", "fst.prime", "gst", "gst.prime", "gst.dbl.prime", "d", "chi2")
    if(stats[1] == "all" & is.haploid(g)) stats <- haploid.stats
    if(stats[1] == "all" & is.diploid(g)) stats <- diploid.stats
    if(is.haploid(g) & !all(stats %in% haploid.stats)) stop("Some 'stats' not valid for haploid data")
    if(is.diploid(g) & !all(stats %in% diploid.stats)) stop("Some 'stats' not valid for diploid data")
    
    list(phist = stat.phist, fst = stat.fst, fst.prime = stat.fst.prime, gst = stat.gst,
         gst.prime = stat.gst.prime.hedrick, gst.dbl.prime = stat.gst.dbl.prime, d = stat.d.jost, chi2 = stat.chi2
    )[stats]
  } else if(is.list(stats)) { # check if list of valid functions
    sapply(stats, function(func) {
      stopifnot.gtype.struct.func(func)
      func
    }, simplify = FALSE)
  } else stop("'stats' must be a character vector or a list of gtype population structure functions")
}