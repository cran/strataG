#' @export gelato.run
#' @aliases gelato
#' @importFrom MASS fitdistr
#' @importFrom parallel mclapply 
#' 
#' @title GELATo - Group ExcLusion and Assignment Test
#' @description Run a GELATo test to evaluate assignment likelihoods of 
#'   groups of samples.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param unknown.strata a character vector listing to assign. Strata must occur in \code{g}.
#' @param nrep number of permutation replicates for Fst distribution.
#' @param min.sample.size minimum number of samples to use to characterize knowns. If 
#'   the known sample size would be smaller than this after drawing an equivalent
#'   number of unknowns for self-assignment, then the comparison is not done.
#' @param num.cores number of CPU cores to use. Value is passed to \code{\link[parallel]{mclapply}}.
#' 
#' @return A list with the following elements:
#' \tabular{ll}{
#'   \code{assign.prob} \tab a data.frame of assignment probabilities.\cr
#'   \code{likelihoods} \tab a list of likelihoods.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references O'Corry-Crowe et. al. XXXX
#' 
#' @seealso \code{\link{gelato.plot}}
#' 
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5)
#'
#' msat.gelato <- gelato.run(msats, "Offshore.South", nrep = 20)
#' gelato.plot(msat.gelato, "Offshore.South")

gelato.run <- function(g, unknown.strata, nrep = 1000, min.sample.size = 5, num.cores = 1) { 
  stopifnot.gtypes(g)
  
  # Check unknown strata
  if(!is.character(unknown.strata) & !is.vector(unknown.strata)) {
    stop("'unknown.strata' must be a character vector")
  }
  all.strata <- attr(g, "strata.names")[-1]
  unknown.strata <- unique(unknown.strata)
  if(!all(unknown.strata %in% all.strata)) {
    stop("Some 'unknown.strata' could not be found in 'g'")
  }
  
  opt <- options(mc.cores = num.cores)
  
  knowns <- sort(setdiff(all.strata, unknown.strata))
  # loop through every unknown strata
  result <- sapply(unknown.strata, function(unknown) {
    unknown.gtypes <- subset(g, strata = unknown)
    unknown.n <- nrow(unknown.gtypes$genotypes)
    
    # loop through each known population and calculate distribution
    #   of Fst and log-likelihood of membership
    unknown.result <- sapply(knowns, function(known) {
      known.gtypes <- subset.gtypes(g, strata = known)
      if((nrow(known.gtypes$genotypes) - unknown.n) >= min.sample.size) {
        fst.dist <- do.call(rbind, mclapply(1:nrep, function(i) {
          # select samples to self assign
          ran.sample <- sample(rownames(known.gtypes$genotypes), unknown.n)
          # extract gtypes of base known strata
          known.to.keep <- setdiff(rownames(known.gtypes$genotypes), ran.sample)
          known.sample <- subset(known.gtypes, ids = known.to.keep)
          # gtypes for observed Fst
          obs.gtypes <- merge(known.sample, unknown.gtypes)
          # gtypes for null Fst
          null.gtypes <- decode(known.gtypes)
          null.gtypes$genotypes[ran.sample, "strata"] <- rep("gelato.unknown", unknown.n)
          gen.mat <- cbind(id = rownames(null.gtypes$genotypes), null.gtypes$genotypes)
          null.gtypes <- gtypes(gen.mat, dna.seq = null.gtypes$sequences)
          c(obs = stat.fst(obs.gtypes)$estimate, null = stat.fst(null.gtypes)$estimate)
        }))
        fst.dist <- fst.dist[apply(fst.dist, 1, function(x) all(!is.na(x))), ]
        
        if(nrow(fst.dist) < 2) {
          NULL
        } else {
          # summarize Fst distribution
          norm.coefs <- fitdistr(fst.dist[, "null"], "normal")$estimate
          log.Lik <- sum(log(dnorm(fst.dist[, "obs"], norm.coefs[1], norm.coefs[2])), na.rm = T)        
          list(fst.dist = fst.dist, 
               log.Lik.smry = c(log.Lik = log.Lik, mean.nreps = log.Lik / length(log.Lik),
                                median = log(dnorm(median(fst.dist[, "obs"], na.rm = T), norm.coefs[1], norm.coefs[2])), 
                                mean = log(dnorm(mean(fst.dist[, "obs"], na.rm = T), norm.coefs[1], norm.coefs[2]))
               ),
               norm.coefs = norm.coefs
          )
        }
      } else NULL
    }, simplify = F)
    
    # calculate median logLikehood of assignment to each known
    log.Lik <- sapply(unknown.result, function(x) {
      if(is.null(x)) NA else x$log.Lik.smry["median"]
    })
    
    print(log.Lik)
    lik <- exp(log.Lik - max(log.Lik, na.rm = T))
    assign.prob <- lik / sum(lik, na.rm = T) 
    names(assign.prob) <- knowns
    
    list(assign.prob = assign.prob, likelihoods = unknown.result)
  }, simplify = F)
  
  assign.prob <- as.data.frame(t(sapply(result, function(x) x$assign.prob)))
  assign.prob$assignment <- apply(assign.prob, 1, function(x) colnames(assign.prob)[which.max(x)])
  result <- lapply(result, function(x) x$likelihoods)
  options(opt)
  list(assign.prob = assign.prob, likelihoods = result)
}