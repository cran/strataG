#' @export low.freq.subs
#' 
#' @title Low Frequency Substitutions
#' @description Check nucleotide sites for low frequency substitutions.
#' 
#' @param dna.seq a list of DNA sequences.
#' @param min.freq minimum frequency of base to be flagged.
#' @param motif.length length of motif around low frequency base to output.
#' 
#' @return data.frame listing id, site number, and motif around low frequency base call.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.haps)
#' low.freq.subs(dolph.haps)

low.freq.subs <- function(dna.seq, min.freq = 3, motif.length = 10) {  
  stopifnot.aligned(dna.seq)
  
  motif.half <- round(motif.length / 2, 0)
  var.sites <- variable.sites(dna.seq)
  has.min.freq <- apply(var.sites$site.freq, 2, function(site.freq) {
    site.freq <- site.freq[site.freq > 0 & site.freq < min.freq]
    length(site.freq) > 0
  })
  sites.w.min.freq <- var.sites$site.freq[, has.min.freq]
  seq.mat <- do.call(rbind, dna.seq)
  sites.to.check <- lapply(colnames(sites.w.min.freq), function(x) {
    position <- as.numeric(x)
    site.freq <- sites.w.min.freq[, x]
    bases <- names(site.freq)[which(site.freq > 0 & site.freq < min.freq)]
    site <- seq.mat[, position]
    id <- names(site)[which(site %in% bases)]
    to.check <- data.frame(id = id, site = rep(position, length(id)))
    to.check$base <- seq.mat[id, position]      
    to.check$motif <- sapply(id, function(i) {
      start.bp <- max(1, position - motif.half)
      end.bp <- min(ncol(seq.mat), position + motif.half)
      paste(seq.mat[i, start.bp:end.bp], collapse = "")
    })
    to.check
  })
  sites.to.check <- do.call(rbind, sites.to.check)  
  sites.to.check[order(sites.to.check$id), ]
}