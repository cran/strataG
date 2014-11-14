#' @rdname decode

decode.loci <- function(g) {
  stopifnot.gtypes(g)
  
  allele.names <- attr(g, "allele.names")
  loci <- g$genotypes[, -1, drop = FALSE]
  decoded <- if(ncol(loci) == 1) {
    hap.mat <- cbind(allele.names[[1]][as.character(loci)])
    rownames(hap.mat) <- rownames(loci)
    hap.mat
  } else {
    do.call(cbind, lapply(seq(1, ncol(loci), 2), function(i) {
      locus.mat <- matrix(as.character(NA), nrow = nrow(loci), ncol = 2)
      k <- (i + 1) / 2
      locus.mat[, 1] <- allele.names[[k]][as.character(loci[, i])]
      locus.mat[, 2] <- allele.names[[k]][as.character(loci[, i + 1])]
      locus.mat
    }))
  }
  rownames(decoded) <- rownames(loci)
  colnames(decoded) <- colnames(loci)
  decoded
}  