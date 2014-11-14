#' @rdname read.fasta

write.fasta <- function(x, file = "sequences.fasta") {
  if(is.gtypes(x, show.warnings = FALSE)) {
    if(!is.null(x$sequences)) {
      x <- decode.sequences(x)
    } else {
      stop("'x' does not contain a list of sequences")
    }
  }
  if(!is.dna.seq(x)) stop("'x' is not a valid list of sequences")
  write.dna(x, file = file, format = "fasta", nbcol = -1, colsep = "", indent = 0, blocksep = 0)
  invisible(file)
}