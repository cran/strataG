#' @rdname gtypes

check.gtypes <- function(g, ploidy = "any") {  
  if(is.null(getOption("check.gtypes"))) options(check.gtypes = TRUE)
  if(!getOption("check.gtypes")) return(g)
  stopifnot.gtypes(g, ploidy)
  
  # remove samples with strata = NA or is NA for all loci
  g <- na.omit.gtypes(g)
  
  # check that at least two strata are present
  if(length(unique(decode.strata(g))) < 2) stop("Fewer than two strata specified")
  
  g
}