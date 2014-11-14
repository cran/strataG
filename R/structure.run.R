#' @rdname structure.run
#' @aliases structure
#' @export structure.run structure.write structure.read
#' @importFrom parallel mclapply
#' 
#' @title STRUCTURE
#' @description Run STRUCTURE to assess group membership of samples.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param k.range vector of values to for \code{maxpop} in multiple runs. If set to \code{NULL}, a single STRUCTURE run is 
#' conducted with \code{maxpops} groups. If specified, do not also specify \code{maxpops}.
#' @param num.k.rep number of replicates for each value in \code{k.range}.
#' @param label label to use for input and output files
#' @param delete.files logical. Delete all files when STRUCTURE is finished?
#' @param num.cores number of CPU cores to use. Value is passed to \code{\link[parallel]{mclapply}}.
#' @param ... arguments to be passed to \code{structure.write}.
#' @param maxpops number of groups.
#' @param burnin number of iterations for MCMC burnin.
#' @param numreps number of MCMC replicates.
#' @param noadmix logical. No admixture?
#' @param freqscorr logical. Correlated frequencies?
#' @param randomize randomize.
#' @param seed set random seed.
#' @param pop.prior a character specifying which population prior model to use, "locprior" or "usepopinfo".
#' @param locpriorinit parameterizes locprior parameter \emph{r} - how informative the populations are. Only used when \code{pop.prior} = "locprior".
#' @param maxlocprior specifies range of locprior parameter \emph{r}. Only used when \code{pop.prior} = "locprior".
#' @param gensback integer defining the number of generations back to test for immigrant ancestry. Only used when \code{pop.prior} = "usepopinfo".
#' @param migrprior numeric between 0 and 1 listing migration prior. Only used when \code{pop.prior} = "usepopinfo".
#' @param pfrompopflagonly logical. update allele frequencies from individuals specified by \code{popflag}. Only used when \code{pop.prior} = "usepopinfo".
#' @param popflag a vector of integers (0, 1) or logicals identifiying whether or not to use strata information. Only used when \code{pop.prior} = "usepopinfo".
#' @param file name of the output file from STRUCTURE.
#' @param pops vector of population labels to be used in place of numbers in STRUCTURE file.
#'    
#' @return
#' \describe{
#'  \item{structure.run}{a list where each element is a list with results from \code{structure.read} and a vector of the filenames used.}
#'  \item{structure.write}{a vector of the filenames used by STRUCTURE.}
#'  \item{structure.read}{a list containing:
#'    \tabular{ll}{
#'      \code{summary} \tab new locus name, which is a combination of loci in group.\cr
#'      \code{q.mat} \tab data.frame of assignment probabilities for each id.\cr
#'      \code{prior.anc} \tab list of prior ancestry estimates for each individual where population priors were used.\cr
#'      \code{files} \tab vector of input and output files used by STRUCTURE.\cr
#'      \code{label} \tab label for the run.\cr
#'    }
#'  }
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references Pritchard, J.K., M. Stephens, P. Donnelly. 2000. Inference of population structure using 
#' multilocus genotype data. Genetics 155:945-959. \url{http://pritchardlab.stanford.edu/structure.html}
#' 
#' @seealso \code{\link{structure.plot}}, \code{\link{structure.evanno}}, \code{\link{clumpp}} 
#' 
#' @examples
#' \dontrun{
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5, description = "msats")
#' 
#' # Run STRUCTURE
#' sr <- structure.run(msats, k.range = 1:6, num.k.rep = 10)
#' 
#' # Calculate Evanno metrics
#' evno <- structure.evanno(sr)
#' print(evno)
#' 
#' # Run CLUMPP to combine runs for K = 2
#' clumpp <- clumpp.run(sr, k = 3)
#' print(clumpp)
#' 
#' # Plot CLUMPP results
#' structure.plot(clumpp)
#' }

structure.run <- function(g, k.range = NULL, num.k.rep = 1, label = NULL, 
                          delete.files = TRUE, num.cores = 1, ...) {
  
  stopifnot.gtypes(g, "diploid")
  
  if(is.null(label)) label <- "structure.run"
  label <- gsub(" ", ".", label)
  
  # setup k and replicate data.frame to cycle through
  if(is.null(k.range)) k.range <- 1:length(unique(g$genotypes[, "strata"]))
  rep.df <- expand.grid(rep = 1:num.k.rep, k = k.range)
  
  unlink(label, recursive = TRUE, force = TRUE)
  dir.create(label)
  if(!file_test("-d", label)) stop(paste("'", label, "' is not a valid folder.", sep = ""))
  label <- file.path(label, label)
  
  rownames(rep.df) <- paste(label, ".k", rep.df$k, ".r", rep.df$rep, sep = "")
  out.files <- mclapply(rownames(rep.df), function(x) {
    sw.out <- structure.write(g, label = x, maxpops = rep.df[x, "k"], ...)
    files <- sw.out$files
    cmd <- paste("structure -m ", files["mainparams"], 
                 " -e ", files["extraparams"], 
                 " -i ", files["data"], 
                 " -o ", files["out"], 
                 sep = ""
    )
    
    err.code <- system(cmd)
    if(err.code == 127) {
      stop("You do not have STRUCTURE installed.")
    } else if(!err.code == 0) {
      stop(paste("Error running STRUCTURE. Error code", err.code, "returned."))
    }
    
    files["out"] <- paste(files["out"], "_f", sep = "")
    result <- structure.read(files["out"], sw.out$pops)
    
    if(file.exists("seed.txt")) file.remove("seed.txt")
    files <- if(delete.files) NULL else files
    
    result <- c(result, list(files = files, label = basename(x)))
    fname <- paste(x, ".ws.rdata", sep = "")
    save(result, file = fname)
    fname
  }, mc.cores = num.cores)
  
  run.result <- lapply(out.files, function(f) {
    result <- NULL
    load(f)
    result
  })
  names(run.result) <- sapply(run.result, function(x) x$label)
  class(run.result) <- c("structure.result", class(run.result))
  
  if(delete.files) unlink(dirname(label), recursive = TRUE, force = TRUE)
  run.result
}