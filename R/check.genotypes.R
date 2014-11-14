#' @export check.genotypes
#' 
#' @title Run All Genotype QA/QC Tests
#' @description Produces by-locus and by-sample summaries, checks for duplicates, 
#' and runs Hardy-Weinberg Jackknife test.
#' 
#' @param x either a .csv file or data.frame with genotype data. First column has IDs, 
#' second column has stratifications, third column to end has genotypes with two columns per locus.
#' @param label label for output folder and prefix for files.
#' @param num.shared either number of loci or percentage of loci two individuals must share to 
#' be considered duplicate individuals.
#' @param exclude.num,min.hwe.samples,alpha parameters for \code{\link[strataG]{jack.hwe}} 
#' and \code{\link[strataG]{jack.influential}}.
#' @param num.cores number of CPU cores to use. Value is passed to \code{\link[parallel]{mclapply}}.
#' 
#' @note Requires that GENEPOP is installed on the system and accessible from the command line.
#' 
#' @return Nothing is returned, but outputs are placed in a folder specified by \code{label}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

check.genotypes <- function(x, label, num.shared = 0.8, exclude.num = 1, 
                      min.hwe.samples = 5, alpha = 0.05, num.cores = 1) {
  
  if(is.character(x)) if(!file.exists(x)) stop("file 'x' can't be found.") else read.gen.data(x)
  if(!file.exists(label)) dir.create(label)
  label <- file.path(label, label)
  
  g <- if(is.gtypes(x, FALSE)) x else gtypes(x)
  stopifnot.gtypes(g, "diploid")
  
  cat(format(Sys.time(), "%H:%M:%S"), "Summarizing loci", "\n")
  smry.by.locus(g, label)
  cat(format(Sys.time(), "%H:%M:%S"), "Summarizing samples", "\n")
  smry.by.sample(g, label)
  cat(format(Sys.time(), "%H:%M:%S"), "Checking for duplicates", "\n")
  dup.genotypes(g, label, num.shared = num.shared, num.cores = num.cores)
  
  cat(format(Sys.time(), "%H:%M:%S"), "Calculating linkage disequilibrium", "\n")
  strata.g <- strata.split(g)
  for(x in names(strata.g)) {
    lnkg <- linkage.genepop(strata.g[[x]], show.output = FALSE)
    file <- paste(label, "linkage.disequilibrium", x, "csv", sep = ".") 
    write.csv(lnkg, file = file, row.names = FALSE)
  }
  lnkg <- linkage.genepop(g, show.output = FALSE)
  file <- paste(label, ".linkage.disequilibrium.all.strata.csv", sep = "")
  write.csv(lnkg, file = file, row.names = FALSE)
  
  cat(format(Sys.time(), "%H:%M:%S"), "HWE jackknife test", "\n")
  all.jack <- jack.hwe(g, exclude.num = exclude.num, min.hwe.samples = min.hwe.samples, show.progress = TRUE)
  infl <- jack.influential(all.jack, alpha = alpha)
  if(!is.null(infl$influential)) {
    write.csv(infl$influential, paste(label, ".influential.ids.csv", sep = ""), row.names = FALSE)
    write.csv(infl$allele.freqs, paste(label, ".influential.allele.freqs.csv", sep = ""), row.names = FALSE)                                            
  
    pdf(file = paste(label, ".jack.influential.plot.pdf", sep = ""))
    to.plot <- plot(infl, basename(label))
    dev.off()
  } else {
    cat("No influential samples found in HWE jackknife test.\n")
  }
  
  cat(format(Sys.time(), "%H:%M:%S"), "Done!\n")
  invisible(NULL)
}
