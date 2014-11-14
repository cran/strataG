#' @export align.folder
#' 
#' @title Align Sequences in a Folder
#' @description Aligns all sequences ending in ".fa*" within a folder
#' 
#' @param input.folder a folder containing sequences.
#' @param output.folder older where alignments will go.
#' @param align.func an alignment function like \code{\link{align.clustal}} or \code{\link{align.mafft}}.
#' @param ... parameters to \code{align.func}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

align.folder <- function(input.folder, output.folder = file.path(input.folder, "Aligned"), align.func, ...) {
  if(!file.exists(output.folder)) dir.create(output.folder) 
  for(fname in dir(input.folder, pattern = ".fa", full.names = TRUE)) {
    dna.seq <- read.dna(fname, format = "fasta", as.character = TRUE, as.matrix = FALSE)
    seqs.aligned <- align.func(dna.seq, ...)
    fname <- file.path(output.folder, paste("aligned", basename(fname)))
    write.dna(seqs.aligned, file = fname, format = "fasta")
  }
}