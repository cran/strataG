#' @export write.snapp.nexus
#' @importFrom ape write.nexus.data
#' 
#' @title Write NEXUS File for SNAPP
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param file the filename the NEXUS file to output.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

write.snapp.nexus <- function(g, file = "snapp.data.nex") {
  stopifnot.gtypes(g, "diploid")
  
  num.alleles <- num.alleles(g)
  biallelic <- names(num.alleles)[num.alleles == 2]
  if(length(biallelic) == 0) {
    warning("No loci are biallelic. No file written.")
    return(NULL)
  }
  
  g <- subset(g, loci = biallelic)
  result <- sapply(attr(g, "locus.names"), function(locus) {
    locus.data <- g$genotypes[, drop(locus.cols(locus, g))]
    genotypes <- apply(locus.data, 1, function(x) {
      switch(paste(sort(x), collapse = "."), '1.1' = 0, '1.2' = 1, '2.2' = 2, NA)
    })
    as.character(as.numeric(factor(genotypes)) - 1)
  })
  
  result[is.na(result)] <- "?"
  result <- lapply(1:nrow(result), function(i) result[i, ])
  strata <- gsub("[ _]", ".", g$genotypes[, "strata"])
  id <- gsub("[ _]", ".", rownames(g$genotypes))
  names(result) <- paste(strata, id, sep = "_")
  
  write.nexus.data(result, file = file)
  
  snapp.file <- scan(file, what = "character", sep = "\n", quiet = TRUE)
  fmt <- grep("FORMAT", snapp.file)
  snapp.file[fmt] <- "  FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS=\"012\" LABELS=LEFT TRANSPOSE=NO INTERLEAVE=NO;"
  write(snapp.file, file = file)
  
  invisible(do.call(rbind, result))
}