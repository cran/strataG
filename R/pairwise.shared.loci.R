#' @export pairwise.shared.loci
#' 
#' @title Calculate Pairwise Shared Loci
#' @description Calculate number of loci shared between pairs of individuals.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param num.cores number of CPU cores to use. Value is passed to \code{\link[parallel]{mclapply}}.
#' 
#' @return matrix summary of pairwise shared loci.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{num.alleles.shared}}

pairwise.shared.loci <- function(g, num.cores = 1) {
  stopifnot.gtypes(g, "diploid")
  opt <- options(mc.cores = num.cores)
  
  id.pairs <- t(combn(rownames(g$genotypes), 2))
  shared <- do.call(rbind, mclapply(1:nrow(id.pairs), function(i) {
    num.alleles.shared(id.pairs[i, 1], id.pairs[i, 2], g)
  }))
  shared.summary <- do.call(rbind, mclapply(1:nrow(shared), function(i) {  
    num.same <- sum(na.omit(shared[i, ]) == 4) 
    num.not.missing <- sum(!is.na(shared[i, ]))
    prop.same <- num.same / num.not.missing
    c(num.same = num.same, num.not.missing = num.not.missing, prop.same = prop.same)
  }))
  shared.summary <- cbind(as.data.frame(id.pairs, stringsAsFactors = FALSE), shared.summary, shared)
  colnames(shared.summary)[1:2] <- c("id1", "id2")
  
  options(opt)
  shared.summary
}