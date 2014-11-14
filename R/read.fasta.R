#' @rdname read.fasta
#' @export read.fasta write.fasta
#' @aliases fasta
#' @importFrom ape read.dna write.dna
#' 
#' @title Read and Write FASTA
#' @description Read and write FASTA formatted files of sequences.
#' 
#' @param file a FASTA-formatted file of sequences.
#' @param x a list of DNA sequences or a haploid \code{\link{gtypes}} object with sequences. 
#' 
#' @return from \code{read.fasta}, a list of DNA sequences.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{as.dna.seq}}

read.fasta <- function(file) {
  dna.seq <- read.dna(file, format = "fasta", as.character = TRUE, as.matrix = FALSE)
  # replace ?s with Ns and convert to lower-case
  lapply(dna.seq, function(x) tolower(gsub("\\?", "n", x)))
}