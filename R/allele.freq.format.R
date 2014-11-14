#' @export allele.freq.format
#' 
#' @title Compiles and Formats Allele Frequencies
#' 
#' @param x data.frame where first column is sample id and second is locus name
#' @param g a \code{gtypes} object.
#' 
#' @return data.frame of original samples, loci, and formatted alleles and frequencies.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#'
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5)
#' 
#' # select 10 random ids and loci to form input data.frame
#' sample.locus.df <- data.frame(
#'   ids = sample(rownames(msats$genotypes), 10, rep = TRUE),
#'   loci = sample(attr(msats, "locus.names"), 10, rep = TRUE)
#' )
#' 
#' allele.freq.format(sample.locus.df, msats)

allele.freq.format <- function(x, g) {    
  if(ncol(x) != 2) stop("'x' must have two columns.")
  stopifnot.gtypes(g)
  
  freqs <- allele.freqs(g)
  g <- decode(g)
  format.freqs <- sapply(1:nrow(x), function(i) {
    id <- as.character(x[i, 1])
    locus <- as.character(x[i, 2])
    locus.num <- which(attr(g, "locus.names") == locus)
    if(length(locus.num) == 0) return(NA)   
    
    if(!id %in% rownames(g$genotypes)) return(NA)
    genotype <- g$genotype[id, locus.cols(locus.num)]
    
    freq.df <- freqs[[locus]]
    this.allele <- freq.df$allele %in% genotype
    freqs <- round(freq.df[this.allele, "prop"], 3)
    freq.1 <- format(freqs[1], nsmall = 3)
    allele.str <- paste(genotype[1], " (", freq.1, ")", sep = "")
    if(!any(is.na(genotype))) {
      if(genotype[1] != genotype[2]) {
        freq.2 <- format(freqs[2], nsmall = 3)
        allele.str <- paste(allele.str, " / ", genotype[2], " (", freq.2, ")", sep = "")
      }
    }
    allele.str
  })
  cbind(x, allele.freqs = format.freqs)
}