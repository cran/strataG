#' @rdname label.haplotypes

label.haplotypes.gtypes <- function(x, ...) {
  stopifnot.gtypes(x, "haploid")
  opt <- options(stringsAsFactors = FALSE) 
  
  # label and reassign
  new.haps <- label.haplotypes(decode.sequences(x), ...)
  
  # reassign haplotypes
  gen.mat <- as.matrix(x)
  gen.mat[, 3] <- new.haps$haps[as.character(gen.mat[, 3])]
  
  options(opt)
  gtypes(gen.mat, id.col = 1, strata.col = 2, locus.col = 3, 
         dna.seq = new.haps$hap.seqs, description = attr(x, "description")
  )
}
