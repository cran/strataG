#' @rdname gtypes

stopifnot.gtypes <- function(g, ploidy = "any") {
  if(is.null(getOption("check.gtypes"))) options(check.gtypes = TRUE)
  if(!getOption("check.gtypes")) invisible(NULL)
  
  # check class
  if(!is.gtypes(g)) stop("'g' is not a valid 'gtypes' class object")
  
  # check ploidy
  if(!ploidy %in% c("any", "haploid", "diploid")) stop("Unknown ploidy type")
  num.cols <- ncol(g$genotypes) - 1
  if(ploidy == "any") if(num.cols != 1 & num.cols %% 2 != 0) stop("'g' is not haploid or diploid")
  if(ploidy == "haploid") if(num.cols != 1) stop("'g' is not haploid")
  if(ploidy == "diploid") if(num.cols %% 2 != 0) stop("'g' is not diploid")
}