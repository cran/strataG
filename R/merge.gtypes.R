#' @export merge.gtypes
#' 
#' @title Merge gtypes 
#' @description Merge two conforming \code{\link{gtypes}} into one.
#' 
#' @usage \method{merge}{gtypes}(x, y, description, ...)
#' 
#' @param x,y two \code{\link{gtypes}} objects.
#' @param description optional description for resulting object. If \code{NULL} then descriptions of both objects are combined.
#' @param ...  other parameters (ignored).
#'  
#' @details \code{x} and \code{y} must have the same ploidy. If haploid, they cannot share samples. If diploid, and shared samples 
#'   cannot also share loci in the two objects.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

merge.gtypes <- function(x, y, description = NULL, ...) {  
  stopifnot.gtypes(x)
  stopifnot.gtypes(y)
  both.haploid <- is.haploid(x) & is.haploid(x)
  both.diploid <- is.diploid(y) & is.diploid(y)
  if(!(both.haploid | both.diploid)) stop("gtypes 'x' and 'y' must have same ploidy")
  
  desc <- if(is.null(description)) paste(attr(x, "description"), "/", attr(y, "description")) else description
  x.df <- as.data.frame(x)
  y.df <- as.data.frame(y)
  if(length(intersect(x.df$id, y.df$id)) != 0) stop("gtypes 'x' and 'y' can't share samples.")
  if(!setequal(colnames(x.df), colnames(y.df))) stop("gtypes 'x' and 'y' must have the same loci.")
  
  if(both.haploid) {
    #x.df$haplotype <- paste(x.df$haplotype, ".x", sep = "")
    x.seq <- decode.sequences(x)
    #if(!is.null(x.seq)) names(x.seq) <- paste(names(x.seq), ".x", sep = "")
    #y.df$haplotype <- paste(y.df$haplotype, ".y", sep = "")
    y.seq <- decode.sequences(y)
    #if(!is.null(y.seq)) names(y.seq) <- paste(names(y.seq), ".y", sep = "")
    gen.mat <- rbind(x.df, y.df)
    dna.seq <- c(x.seq, y.seq)
    gtypes(gen.mat, dna.seq = dna.seq, description = desc)    
  } else {
    gen.mat <- rbind(x.df, y.df)
    locus.data <- gen.mat[, -c(1:2)]
    locus.data <- locus.data[, order(colnames(locus.data))]
    gen.mat <- cbind(gen.mat[, 1:2], locus.data)
    gtypes(gen.mat, description = desc)
  }
}