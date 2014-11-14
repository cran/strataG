#' @export relabel.sequences
#' 
#' @title Relabel Sequences
#' @description Relabel sequences in a \code{gtypes} object by changing names or adding to names.
#' 
#' @param g \code{\link{gtypes}} object.
#' @param new.labels character vector giving new labels for sequences. Must either be of length 1 or the 
#'   number of sequences in \code{g}. If of length 1 and \code{add.labels} = FALSE, sequences will be
#'   sequentially numbered.
#' @param add.labels logical. Should \code{new.labels} be added to the existing sequence labels?
#' @param sep character separator between new and old labels (if \code{add.labels} = FALSE), or new 
#'   label and sequential numbers.
#' 
#' @return a \code{gtypes} object with sequences relabelled.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' @seealso \code{\link{as.dna.seq}}

relabel.sequences <- function(g, new.labels, add.labels = TRUE, sep = ".") {
  stopifnot.gtypes(g, "haploid")
  
  if(is.null(g$sequences)) stop("'g' does not contain sequences.")
  new.labels <- as.character(new.labels)
  if(!length(new.labels) %in% c(1, length(g$sequences))) stop("'new.labels' must be of length 1 or equal to the number of sequences.")
  if(length(new.labels) != length(unique(new.labels))) stop("All 'new.labels' must be unique.")
  
  desc <- attr(g, "description")
  g <- decode(g)
  
  old.labels <- names(g$sequences)
  new.labels <- if(add.labels) {
    paste(new.labels, old.labels, sep = sep)
  } else if(length(new.labels) == length(g$sequences)) {
    new.labels
  } else {
    paste(new.labels, 1:length(old.labels), sep = sep)
  }
  
  names(g$sequences) <- new.labels
  names(new.labels) <- old.labels
  gen.mat <- g$genotypes
  gen.mat[, 2] <- new.labels[gen.mat[, 2]]
  gen.mat <- cbind(id = rownames(gen.mat), gen.mat)
  
  gtypes(gen.mat, dna.seq = g$sequences, description = desc)
}