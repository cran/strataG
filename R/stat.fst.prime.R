#' @export stat.fst.prime
#' 
#' @title F'st analysis of population structure
#' 
#' @param g a \code{gtypes} object
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5, description = "msats")
#' 
#' stat.fst.prime(msats)

stat.fst.prime <- function(g) {  
  g <- check.gtypes(g, "diploid") 
  
  fst.max.g <- decode(g)
  fst.max.g$genotypes[, -1] <- t(sapply(1:nrow(fst.max.g$genotypes), function(i) {
    na.locs <- is.na(fst.max.g$genotypes[i, -1])
    new.locus <- paste(fst.max.g$genotypes[i, 1], fst.max.g$genotypes[i, -1], sep = ".")
    new.locus[na.locs] <- NA
    new.locus
  }))
  fst.max.g <- gtypes(cbind(id = rownames(fst.max.g$genotypes), fst.max.g$genotypes))
  
  fst.max <- stat.fst(fst.max.g)
  fst <- stat.fst(g)
  fst$estimate <- fst$estimate / fst.max$estimate
  fst$stat.name <- "F'st"
  fst
}
class(stat.fst.prime) <- c(class(stat.fst.prime), "gtype.struct.func")