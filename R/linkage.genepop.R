#' @export linkage.genepop
#' 
#' @title Linkage Disequlibrium 
#' @description Calculate linkage disequilibrium p-values using GENEPOP.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param show.output logical. Show GENEPOP output on console?
#' @param delete.files logical. Delete GENEPOP input and output files when done?
#' @param label character string to use to label GENEPOP input and output files.
#' @param ... other arguments to be passed to \code{\link{genepop.run}}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

linkage.genepop <- function(g, show.output = TRUE, delete.files = TRUE, label = "linkage.genepop", ...) {
  stopifnot.gtypes(g, "diploid")
  warn <- NULL
  opt <- options(warn = -1)
  
  # Run Genepop
  g <- restratify(g, rep("1", nrow(g$genotypes)))
  output <- genepop.run(g, output.ext = ".DIS", show.output = show.output, label = label, other.settings = "MenuOptions=2.1")
  if(!is.list(output)) {
    options(opt)
    return(NA)
  }
  result <- scan(output$files["output.fname"], what = "character", quiet = TRUE)
  loc.names <- output$locus.names
  
  # Create empty matrix
  numrows <- ((length(loc.names) ^ 2) - length(loc.names)) / 2
  result.mat <- matrix(as.character(NA), numrows, 5)
  
  # Find starting points
  loc <- grep("Switches", result, value = F) + 7
  first.col <- grep(names(loc.names)[1], result)[1]
  num.skip <- first.col - loc
  row.mask <- c(rep(F, num.skip), rep(T, 5))
  
  # Read matrix
  for(r in 1:numrows) {
    result.mat[r, ] <- result[loc:(loc + num.skip + 5)][row.mask]
    loc <- loc + num.skip + 5
  }
  
  # Convert to data.frame and format columns
  result.df <- data.frame(result.mat, stringsAsFactors = FALSE)
  colnames(result.df) <- c("Locus.1", "Locus.2", "p.value", "std.err", "switches")
  result.df$p.value <- as.double(result.df$p.value) #no longer part of R because no longer needed
  result.df$std.err <- as.double(result.df$std.err) # changed from "as.real()"
  result.df$switches <- as.integer(result.df$switches)
  result.df$Locus.1 <- loc.names[result.df$Locus.1]
  result.df$Locus.2 <- loc.names[result.df$Locus.2]    
  
  if(delete.files) for(f in output$files) if(file.exists(f)) file.remove(f) 
  options(opt)
  result.df
}