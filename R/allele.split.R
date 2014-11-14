#' @export allele.split
#' 
#' @title Split Alleles For Diploid Data
#' @description Split loci from one column to two columns in a matrix of diploid data.

#' @param x a matrix or data.frame containing diploid data. Every column represents one locus with two alleles.
#' @note Alleles should be of equal length (e.g., 145095 = 145 and 095, or AG is A and G).
#' 
#' @return matrix with alleles for each locus in one column split into separate columns.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

allele.split <- function(x) {
  opt <- options(warn = -1)
  locus.names <- if(is.null(colnames(x))) paste("Locus", 1:ncol(x), sep = "") else colnames(x)
  x <- do.call(rbind, lapply(1:nrow(x), function(r) as.character(x[r, ])))
 
  split.alleles <- apply(x, 1, function(genotype) {
    genotype.split <- lapply(genotype, function(locus) {
      locus <- sub(" ", "", locus)
      if (is.na(locus)) return(c(NA, NA))
      end <- nchar(locus)
      half <- end / 2
      
      a1 <- substr(locus, 1, half)
      a1 <- ifelse(is.na(as.numeric(a1)), a1, ifelse(as.numeric(a1) == 0, NA, a1))
      a1 <- ifelse((a1 == "NA") | (a1 == ""), NA, a1)
      
      a2 <- substr(locus, half + 1, end)
      a2 <- ifelse(is.na(as.numeric(a2)), a2, ifelse(as.numeric(a2) == 0, NA, a2))
      a2 <- ifelse((a2 == "NA") | (a2 == ""), NA, a2)
      
      c(a1, a2)
    })
    unlist(genotype.split)
  })
  split.alleles <- t(split.alleles)
  colnames(split.alleles) <- paste(rep(locus.names, each = 2), c(1, 2), sep = ".")
  rownames(split.alleles) <- NULL
  options(opt)
  split.alleles
}