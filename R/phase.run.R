#' @rdname phase.run
#' @aliases phase
#' @export phase.run phase.write phase.read.sample phase.read.pair phase.posterior phase.filter
#' @importFrom parallel mclapply 
#' 
#' @title PHASE
#' @description Run PHASE to estimate the phase of loci in diploid data.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param loci vector of names of locus in 'gtypes' that are to be phased.
#' @param positions position along chromosome of each locus.
#' @param type type of each locus.
#' @param num.iter number of PHASE MCMC iterations.
#' @param thinning number of PHASE MCMC iterations to thin by.
#' @param burnin number of PHASE MCMC iterations for burnin.
#' @param model PHASE model type.
#' @param ran.seed PHASE random number seed.
#' @param final.run.factor optional.
#' @param save.posterior logical. Save posterior sample in output list?
#' @param in.file name to use for PHASE input file.
#' @param out.file name to use for PHASE output files.
#' @param delete.files logical. Delete PHASE input and output files when done?
#' @param ph.res result from \code{phase.run}.
#' @param thresh minimum probability for a genotype to be selected (0.5 - 1).
#' @param keep.missing logical. T = keep missing data from original data set. F = Use estimated genotypes from PHASE.
#' @param num.cores number of CPU cores to use. Value is passed to \code{\link[parallel]{mclapply}}.
#'  
#' @note Assumes that the the command line version of PHASE is properly installed and available on the command line,
#'   so it is executable from any directory. On PC's, this requires having it in a folder in 
#'   the PATH environmental variable. On Macs, the executable should be installed in a folder 
#'   like \code{/usr/local/bin} 
#' 
#' @details
#' \tabular{ll}{
#'   \code{phase.run} \tab runs PHASE assuming that the executable is installed properly and available on the command line.\cr
#'   \code{phase.write} \tab writes a PHASE formatted file.\cr
#'   \code{phase.read.pair} \tab reads the '_pair' output file.\cr
#'   \code{phase.read.sample} \tab reads the '_sample' output file.\cr
#'   \code{phase.filter} \tab filters the result from \code{phase.run} to extract one genotype for each sample.\cr
#'   \code{phase.posterior} \tab create a data.frame all genotypes for each posterior sample.\cr
#' }
#'  
#' @return
#' \describe{
#'  \item{phase.run}{a list containing:
#'    \tabular{ll}{
#'      \code{locus.name} \tab new locus name, which is a combination of loci in group.\cr
#'      \code{gtype.probs} \tab a data.frame listing the estimated genotype for every sample along with probability.\cr
#'      \code{orig.gtypes} \tab the original gtypes object for the composite loci.\cr
#'      \code{posterior} \tab a list of \code{num.iter} data.frames representing posterior sample of genotypes for each sample.\cr
#'    }}
#'  \item{phase.write}{a list with the input filename and the \code{\link{gtypes}} object used.}
#'  \item{phase.read.pair}{a data.frame of genotype probabilities.}
#'  \item{phase.read.sample}{a list of data.frames representing the posterior sample of genotypes for one set of loci for each sample.}
#'  \item{phase.filter}{a matrix of genotypes for each sample.}
#'  \item{phase.posterior}{a list of data.frames representing the posterior sample of all genotypes for each sample.}
#' }
#' 
#' @references Stephens, M., and Donnelly, P. (2003). A comparison of Bayesian methods for haplotype reconstruction from 
#'   population genotype data. American Journal of Human Genetics 73:1162-1169.
#'   Available at: \url{http://stephenslab.uchicago.edu/software.html#phase}
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples \dontrun{
#' data(bowhead.snps)
#' data(bowhead.snp.position)
#' snps <- gtypes(bowhead.snps, description = "Bowhead SNP example")
#' summary(snps)
#' 
#' # Run PHASE on all data
#' phase.results <- phase.run(snps, bowhead.snp.position, num.iter = 100, save.posterior = FALSE)
#' 
#' # Filter phase results
#' filtered.results <- phase.filter(phase.results, thresh = 0.5)
#' 
#' # Convert phased genotypes to gtypes
#' ids <- rownames(filtered.results)
#' strata <- bowhead.snps$Stock[match(ids, bowhead.snps$LABID)]
#' filtered.results <- cbind(id = ids, strata = strata, filtered.results)
#' phased.snps <- gtypes(filtered.results, description = "Bowhead phased SNPs")
#' summary(phased.snps)
#' }

phase.run <- function(g, loci, positions = NULL, type = NULL,
  num.iter = 100000, thinning = 100, burnin = 100000, model = "new", 
  ran.seed = NULL, final.run.factor = NULL, save.posterior = FALSE, 
  in.file = "phase_in", out.file = "phase_out", delete.files = TRUE, num.cores = 1) { 
  
  stopifnot.gtypes(g, "diploid")
  
  if(!is.data.frame(loci)) {
    if(!(is.character(loci) & is.vector(loci))) stop("'loci' must be a data.frame or character vector")
    if(is.null(positions)) positions <- rep(1, length(locus.names))
    if(length(positions) != length(loci)) stop("'positions' must be same length as 'loci'")
    loci <- data.frame(locus = loci, position = positions, group = 1)
  }
  loci$group <- as.character(loci$group)
  loci$position <- as.numeric(loci$position)
  
  if(is.null(type)) type <- rep("S", length(unique(loci$group)))
  if(length(type) != length(unique(loci$group))) stop("'type' must be same length as number of locus groups")
  names(type) <- unique(loci$group)
  
  opt <- options(mc.cores = num.cores)
  result <- mclapply(unique(loci$group), function(grp) {
    ran.lets <- paste(sample(c(0:9, letters), 10, replace = TRUE), collapse = "")
    in.file <- paste("phase_in_", ran.lets, sep = "")
    out.file <- paste("phase_out_", ran.lets, sep = "")
    
    # Write input file
    group.df <- loci[loci$group == grp, ]  
    locus.type <- rep(type[grp], nrow(group.df))
    in.file.data <- phase.write(g, loci = group.df$locus, positions = group.df$position, type = locus.type, in.file)

    # Set parameters
    M.opt <- switch(model, new = "-MR", old = "-MS", hybrid = "-MQ", "")           
    S.opt <- ifelse(is.null(ran.seed), "", paste("-S", ran.seed, sep = ""))
    X.opt <- ifelse(is.null(final.run.factor), "", paste("-X", final.run.factor, sep = ""))
    s.opt <- ifelse(save.posterior, "-s", "")
    in.file.opt <- paste("\"", in.file, "\"", sep = "")
    out.file.opt <- paste("\"", out.file, "\"", sep = "")
    iter.params <- paste(trunc(num.iter), trunc(thinning), trunc(burnin))
    phase.cmd <- paste("PHASE", M.opt, S.opt, X.opt, s.opt, in.file.opt, out.file.opt, iter.params)
    
    # Run Phase
    err.code <- system(phase.cmd)  
    if(err.code == 127) {
      stop("You do not have PHASE installed.") 
    } else if(!err.code == 0) {
      stop(paste("Error running PHASE. Error code", err.code, "returned."))
      cat("\n")
    }
    
    # Read output
    gtype.probs <- phase.read.pair(paste(out.file, "_pairs", sep = ""))
    if(is.null(gtype.probs)) {
      alleles <- rep(NA, nrow(g$genotypes))
      gtype.probs <- data.frame(id = rownames(g$genotypes), a1 = alleles, a2 = alleles, pr = rep(1, nrow(g$genotypes)))
    }
    new.locus.name <- paste(group.df$locus, collapse = "_")
    alleles <- paste(new.locus.name, 1:2, sep = ".")
    colnames(gtype.probs)[1:3] <- c("id", alleles) 
    rownames(gtype.probs) <- NULL
    
    locus.result <- list(locus.name = new.locus.name, gtype.probs = gtype.probs, orig.gtypes = in.file.data$gtypes)  
    
    if(save.posterior) {
      locus.result$posterior <- phase.read.sample(paste(out.file, "_sample", sep = ""), paste(locus.type, collapse = ""))
      for(i in 1:length(locus.result$posterior)) colnames(locus.result$posterior[[i]]) <- c("id", alleles)
    }
    
    if(delete.files) file.remove(c(dir(pattern = in.file), dir(pattern = out.file)))
    
    locus.result  
  })
  options(opt)
  names(result) <- lapply(result, function(x) x$locus.name)
  class(result) <- c("phase.result", class(result))
  result
}