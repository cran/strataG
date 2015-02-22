#' @rdname gtypes

is.gtypes <- function(x, show.warnings = TRUE) {  
  if(is.null(getOption("check.gtypes"))) options(check.gtypes = TRUE)
  if(!getOption("check.gtypes")) return(TRUE)
  
  if(!(is.list(x) & inherits(x, "gtypes"))) {
    if(show.warnings) warning("'x' is of the wrong class.")
    return(FALSE)
  }

  w <- ""
  if(!is.null(x$genotypes)) {
    genotypes.check <- is.numeric(x$genotypes) & is.matrix(x$genotypes)
    if(genotypes.check) {
      strata.col <- which(colnames(x$genotypes) == "strata")
      strata.check <- if(length(strata.col) > 0) strata.col == 1 else FALSE
      loci.check <- ncol(x$genotypes) == 2 | (ncol(x$genotypes[, -1, drop = FALSE]) %% 2) == 0
      if(strata.check & loci.check) {
        locus.names <- attr(x, "locus.names")
        has.locus.names <- !is.null(locus.names) 
        if(has.locus.names) {
          locus.name.check <- length(locus.names) == num.loci(x$genotypes[, -1, drop = FALSE])
          haplotype.name.check <- if(!is.null(x$sequences)) {
            all(na.omit(unique(as.character(x$genotypes[, 2]))) %in% names(x$sequences))
          } else TRUE
          if(locus.name.check & haplotype.name.check) {
            return(TRUE)
          } else {
            if(!haplotype.name.check) w <- c(w, "There are haplotypes in 'x' that are missing sequences.")
            if(!locus.name.check) w <- c(w, "The length of the 'locus.names' attributes does not match the number of loci in 'x'.")
          }
        } else w <- c(w, "There is no 'locus.names' attribute in 'x'.")
      } else {
        if(!strata.check) w <- c(w, "There is either no column named 'strata' in the 'genotypes' element of 'x', or it is not the first column.")
        if(!loci.check) w <- c(w, "The number of genotype data columns in 'x' is not 1 or an even number.")
      }
    } else w <- c(w, "The 'genotypes' element in 'x' is not a numeric matrix.")
  } else w <- c(w, "There is no 'genotypes' element in 'x'.")

  if(show.warnings) warning(paste(w, sep = "\n"))      
  FALSE     
}