#' @export check.seqs
#' 
#' @title Run All Sequence QA/QC Tests
#' @description Produces haplotypic summaries, checks for low-frequency substitutions, 
#'   and unusually distant sequences.
#' 
#' @param x either a .fasta file or list of sequences.
#' @param label label for output folder and prefix for files.
#' @param min.freq minimum frequency of base to be flagged.
#' @param motif.length length of motif around low frequency base to output.
#' @param model a character string specifying the evolutionary model to be used. Passed to 
#'   \code{\link[ape]{dist.dna}}. 
#' @param pairwise.deletion a logical indicating whether to delete the sites with missing data 
#'   in a pairwise way. Passed to \code{\link[ape]{dist.dna}}. 
#' 
#' @return Nothing is returned, but outputs are placed in a folder specified by \code{label}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

check.seqs <- function(x, label, min.freq = 3, motif.length = 10, model = "raw", pairwise.deletion = FALSE) {  
  if(is.character(x)) if(!file.exists(x)) stop("file 'x' can't be found.") else x <- read.fasta(x)
  if(!is.dna.seq(x)) stop("'x' is not a set of DNA sequences.")
  if(!file.exists(label)) dir.create(label)
  label <- file.path(label, label)
  
  dna.seq <- if(is.aligned.seq(x)) x else align.mafft(x)
  
  cat(format(Sys.time(), "%H:%M:%S"), "Checking haplotypes", "\n")
  hap.check <- haplotype.likelihoods(dna.seq, model = model, pairwise.deletion = pairwise.deletion)
  hap.check <- cbind(hap.check)
  colnames(hap.check) <- "hap.delta.lnL"
  file = paste(label, ".hap.likelihoods.csv", sep = "")
  write.csv(cbind(hap.check), file = file)
  
  cat(format(Sys.time(), "%H:%M:%S"), "Checking sites", "\n")
  site.check <- low.freq.subs(dna.seq, min.freq = min.freq, motif.length = motif.length)
  if(nrow(site.check) != 0) {
    file = paste(label, ".site.check.csv", sep = "")
    write.csv(site.check, file = file, row.names = FALSE)
  }
    
  cat(format(Sys.time(), "%H:%M:%S"), "Done!\n")
  invisible(NULL)
}