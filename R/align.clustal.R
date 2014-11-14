#' @export align.clustal
#' @importFrom ape write.dna
#' @importFrom ape read.dna
#' @importFrom ape read.nexus.data
#' 
#' @title CLUSTAL Alignment
#' @description Align sequences using CLUSTAL
#' 
#' @param x a \code{\link{gtypes}} object with aligned sequences, or a list of aligned DNA sequences.
#' @param gapopen penalty for gap opening.
#' @param gapext penalty for gap extension.
#' @param trans.weight transition weighting.
#' @param input.name name of input file (without extension).
#' @param output.format format of output file ("CLUSTAL", "PHYLIP", "NEXUS").
#' @param output.name name of output file (without extension).
#' @param delete.output delete the output file when finished?
#' 
#' @note Formats and executes a call to the executable \emph{clustalw2}, assuming that it is installed
#'   on the system and available at the command line. 
#' 
#' @return list of aligned sequences
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references CLUSTALW2 is available at: \url{http://www.clustal.org/clustal2}

align.clustal <- function(x, gapopen = 5, gapext = 15, trans.weight = 0.5,
                          input.name = "to_align", output.format = "CLUSTAL", output.name = "aligned_seqs",
                          delete.output = TRUE) {
  if(is.gtypes(x, FALSE)) x <- decode.sequences(x)
  if(!is.dna.seq(x, FALSE)) stop("'x' must contain valid sequences")
  
  input.fasta <- paste(input.name, ".fasta", sep = "")
  infile <- paste("-INFILE=", input.fasta, sep = "")
  g.open <- paste("-GAPOPEN=", gapopen, sep = "")
  g.ext <- paste("-GAPEXT=", gapext, sep = "")
  tw <- paste("-TRANSWEIGHT=", trans.weight, sep = "")
  output <- paste("-OUTPUT=", output.format, sep = "")
  out.ext <- switch(output.format,
                    CLUSTAL = ".aln",
                    PHYLIP = ".phy",
                    NEXUS = ".nex"
  )
  output.name <- paste(output.name, out.ext, sep = "")
  outfile <- paste("-OUTFILE=", output.name, sep = "")
  clustal.call <- paste("clustalw2 -ALIGN -TYPE=DNA -OUTORDER=INPUT",
                        infile, g.open, g.ext, tw, output, outfile
  )
  
  write.dna(x, file = input.fasta, format = "fasta")
  err.code <- system(clustal.call, intern = F)
  if(!err.code == 0) return(NA)
  aligned <- switch(output.format,
                    CLUSTAL = read.dna(output.name, format = "clustal", as.character = TRUE, as.matrix = FALSE),
                    PHYLIP = read.dna(output.name, format = "interleaved", as.character = TRUE, as.matrix = FALSE),
                    NEXUS = read.nexus.data(output.name)
  )
  
  file.remove(input.fasta, paste(input.name, ".dnd", sep = ""))
  if(delete.output) file.remove(output.name)
  lapply(as.dna.seq(aligned), tolower)
}