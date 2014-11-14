#' @export write.gtypes
#' 
#' @title write.gtypes
#' @description Write a \code{\link{gtypes}} object to file(s).
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param label label for filename(s). Default is the gtypes description if present.
#' @param folder folder where file(s) should be written to. If NULL, files are written to current working directory.
#' @param as.frequency logical indicating if haploid data should be output as frequency tables.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

write.gtypes <- function(g, label = NULL, folder = NULL, as.frequency = FALSE) {
  stopifnot.gtypes(g)
  
  desc <- attr(g, "description")
  label <- if(!is.null(label)) label else if(!is.null(desc)) desc else "strataG.gtypes"
  fname <- paste(label, ".csv", sep = "")
  if(!is.null(folder)) fname <- file.path(folder, fname)
  csv.file <- if(is.haploid(g) & as.frequency) {
    x <- as.frequency(g) 
    x <- data.frame(haplotype = rownames(x), cbind(x))
    rownames(x) <- NULL
    x
  } else as.matrix(g)
  write.csv(csv.file, file = fname, row.names = FALSE)
  if(!is.null(g$sequences)) {
    fname <- paste(label, ".fasta", sep = "")
    if(!is.null(folder)) fname <- file.path(folder, fname)
    write.fasta(decode.sequences(g), fname)
  }
}