#' @rdname structure.run

structure.write <- function(g, label = NULL, maxpops = length(unique(g$genotypes[, "strata"])), 
  burnin = 1000, numreps = 1000, noadmix = TRUE, freqscorr = FALSE, 
  randomize = TRUE, seed = 0, pop.prior = NULL,
  locpriorinit = 1, maxlocprior = 20, gensback = 2, migrprior = 0.05,
  pfrompopflagonly = TRUE, popflag = NULL, ...) {
  
  stopifnot.gtypes(g, "diploid")
  g <- decode(g)
  
  # check parameters
  if(!is.null(pop.prior)) {
    if(!pop.prior %in% c("locprior", "usepopinfo")) stop("'pop.prior' must be 'locprior' or 'usepopinfo'.")
  }
  if(is.null(popflag)) popflag <- rep(1, nrow(g$genotypes))
  popflag <- as.numeric(popflag)
  if(length(popflag) != nrow(g$genotypes)) stop("'popflag' should be the same length as the number of individuals in 'g'.")
  if(!all(popflag %in% c(0, 1))) stop("all values in 'popflag' must be 0 or 1.")
  
  in.file <- ifelse(is.null(label), "data", paste(label, "data", sep = "_"))
  out.file <- ifelse(is.null(label), "out", paste(label, "out", sep = "_"))
  main.file <- ifelse(is.null(label), "mainparams", paste(label, "mainparams", sep = "_"))
  extra.file <- ifelse(is.null(label), "extraparams", paste(label, "extraparams", sep = "_"))
  
  # write data
  write(paste(attr(g, "locus.names"), collapse = " "), file = in.file)
  pop.fac <- factor(g$genotypes[, "strata"])
  popdata <- as.numeric(pop.fac)
  
  locus.data <- g$genotypes[, -1, drop = FALSE]
  locus.data[is.na(locus.data)] <- -9
  for(i in 1:nrow(g$genotypes)) {
    write(paste(rownames(g$genotypes)[i], popdata[i], popflag[i], paste(locus.data[i, ], collapse = " ")), 
          file = in.file, append = TRUE
    )
  }
  
  # write mainparams
  main.params <- c(
    paste("MAXPOPS", as.integer(maxpops)),
    paste("BURNIN", as.integer(burnin)),
    paste("NUMREPS", as.integer(numreps)),
    paste("INFILE", in.file),
    paste("OUTFILE", out.file),
    paste("NUMINDS", nrow(g$genotypes)),
    paste("NUMLOCI", length(attr(g, "locus.names"))),
    "MISSING -9",
    "ONEROWPERIND 1",
    "LABEL 1",
    "POPDATA 1",
    "POPFLAG 1",
    "LOCDATA 0",
    "PHENOTYPE 0",
    "EXTRACOLS 0",
    "MARKERNAMES 1"
  )
  main.params <- paste("#define", main.params)
  write(main.params, file = main.file)
  
  # write extraparams
  extra.params <- c(
    paste("NOADMIX", as.integer(noadmix)),
    paste("FREQSCORR", as.integer(freqscorr)),
    "INFERALPHA 1",
    "ALPHA 1.0",
    "FPRIORMEAN 0.01",
    "FPRIORSD 0.05",
    "LAMBDA 1.0",
    "UNIFPRIORALPHA 1", 
    "ALPHAMAX 20.0",
    "ALPHAPRIORA 0.05",
    "ALPHAPRIORB 0.001",
    "COMPUTEPROB 1",
    paste("ADMBURNIN", max(0, as.integer(burnin / 2))),
    "ALPHAPROPSD 0.025",
    "STARTATPOPINFO 0",
    paste("RANDOMIZE", as.integer(randomize)),
    paste("SEED", as.integer(seed)),
    "METROFREQ 10",
    "REPORTHITRATE 0" 
  )
  
  if(!is.null(pop.prior)) {
    pop.prior <- tolower(pop.prior)
    prior.params <- if(pop.prior == "locprior") {
      c("LOCPRIOR 1",
        "LOCISPOP 1",
        paste("LOCPRIORINIT", locpriorinit),
        paste("MAXLOCPRIOR", maxlocprior)
      )        
    } else if(pop.prior == "usepopinfo") {
      c("USEPOPINFO 1",
        paste("GENSBACK", trunc(gensback)),
        paste("MIGRPRIOR", migrprior),
        paste("PFROMPOPFLAGONLY", as.integer(pfrompopflagonly))
      )
    }
    extra.params <- c(extra.params, prior.params)
  }

  extra.params <- extra.params[!is.na(extra.params)]
  extra.params <- paste("#define", extra.params)
  write(extra.params, file = extra.file)
  
  invisible(list(
    files = c(data = in.file, mainparams = main.file, extraparams = extra.file, out = out.file),
    pops = levels(pop.fac)
  ))
}
