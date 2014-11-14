#' @rdname fastsimcoal.run
#' @export fastsimcoal.run fastsimcoal.write
#' @aliases fastsimcoal
#' 
#' @title Run FASTSIMCOAL
#' @description Run FASTSIMCOAL to generate a list of simulated gtypes.
#' 
#' @param num.pops number of populations.
#' @param Ne effective population size.
#' @param sample.size number of samples to take.
#' @param sample.time time to draw samples.
#' @param growth.rate growth rate of populations.
#' @param mig.mat migration matrix.
#' @param hist.ev historical events.
#' @param num.chrom number of chromosomes.
#' @param data.type type of data.
#' @param locus.params locus parameters.
#' @param label character string to label files with.
#' @param num.sims number of simulations to run.
#' @param inf.site.model logical. Infinite site model?
#' @param quiet logical. Run quietly?
#' @param delete.files logical. Delete files when done?
#' 
#' @note Assumes that the program \code{fastsimcoal} is properly installed and available on the command line.
#'   On PC's, this requires having it in a folder in 
#'   the PATH environmental variable. On Macs, the executable should be installed in a folder 
#'   like \code{/usr/local/bin} 
#' 
#' @return A list of \code{\link{gtypes}} objects for each simulated dataset.
#' 
#' @references Excoffier, L. and Foll, M (2011) fastsimcoal: a continuous-time coalescent 
#'   simulator of genomic diversity under arbitrarily complex evolutionary scenarios 
#'   Bioinformatics 27: 1332-1334.  \url{http://cmpg.unibe.ch/software/fastsimcoal2/}
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

fastsimcoal.run <- function(num.pops, Ne, sample.size, sample.time = rep(0, num.pops),
                            growth.rate = rep(0, num.pops), mig.mat = NULL, 
                            hist.ev = NULL, num.chrom = 1, data.type = NULL,
                            locus.params = NULL, label = "fastsimcoal",
                            num.sims = 1, inf.site.model = F, quiet = T, delete.files = T) {
  
  hist.ev <- if(is.list(hist.ev)) do.call(rbind, hist.ev) else rbind(hist.ev)
  locus.params <- if(is.list(locus.params)) do.call(rbind, locus.params) else rbind(locus.params)
  if(nrow(locus.params) == 1 & num.chrom > 1) locus.params <- do.call(rbind, lapply(1:num.chrom, function(i) locus.params[1, ]))
  
  # Write input file
  file <- fastsimcoal.write(num.pops = num.pops, Ne = Ne, sample.size = sample.size, 
    sample.time = sample.time, growth.rate = growth.rate, mig.mat = mig.mat, 
    hist.ev = hist.ev, num.chrom = num.chrom, data.type = data.type, 
    locus.params = locus.params, label = label
  )
  
  # Run fastsimcoal
  if(file.exists(label)) for(f in dir(label, full.names = T)) file.remove(f)
  cmd <- paste("fastsimcoal -i", file, "-n", num.sims, 
    ifelse(inf.site.model, "-I", ""), ifelse(quiet, "-q", ""), sep = " "
  )
  err <- system(cmd, intern = F)
  
  # Read arlequin output
  arl.files <- dir(label, pattern = ".arp", full.names = T)
  result <- if(err != 0 | length(arl.files) == 0) {
    if(err != 0) stop(paste("fastsimcoal returned error code", err)) else stop("fastsimcoal did not generate output")
  } else {
    arl.gtypes <- lapply(arl.files, function(file) {
      f <- readLines(file)
      
      # get start and end points of data blocks
      start <- grep("SampleData=", f) + 1
      end <- which(f == "}") - 2
      pos <- cbind(start, end)
      # compile data for each population
      pop.data <- do.call(rbind, lapply(1:nrow(pos), function(i) {
        f.line <- f[pos[i, 1]:pos[i, 2]]
        f.line <- gsub("[[:space:]]+", "--", f.line)
        data.mat <- do.call(rbind, strsplit(f.line, "--"))[, -2]
        data.mat <- cbind(rep(paste("Sample", i), nrow(data.mat)), data.mat)
      }))
      
      # get data type
      data.type <- f[grep("DataType=", f)]
      data.type <- gsub("\tDataType=", "", data.type)
      is.haploid <- switch(data.type, DNA = T, MICROSAT = F, STANDARD = F)
      if(is.haploid) {      
        # replace sequence with all A's if there are no variable sites
        n.loc <- locus.params[1, 1]
        if(pop.data[1, 3] == "?") {
          full.seq <- paste(rep("A", n.loc), collapse = "")
          pop.data[, 3] <- rep(full.seq, nrow(pop.data))
        } else { # otherwise add A's to pad out to full sequence length
          partial.seq <- paste(rep("A", n.loc - nchar(pop.data[1, 3])), collapse = "")
          pop.data[, 3] <- sapply(pop.data[, 3], function(x) paste(x, partial.seq, sep = "", collapse = ""))
        }
        dna.seq <- strsplit(pop.data[, 3], "")
        names(dna.seq) <- pop.data[, 2]
        g <- gtypes(dna.seq, strata = pop.data[, 1], description = file)
        label.haplotypes(g, "Hap.")
      } else {
        # compile diploid data
        n.loc <- ncol(pop.data) - 2
        pop.data <- do.call(rbind, lapply(seq(1, nrow(pop.data), 2), function(i) {
          ind <- pop.data[c(i, i + 1), ]
          locus.data <- as.vector(ind[, -(1:2)])
          c(ind[1, 1], paste(ind[, 2], collapse = "/"), locus.data)
        }))
        pop.data <- data.frame(pop.data[, 2], pop.data[, 1], pop.data[, -(1:2)])
        gtypes(pop.data, description = file)
      }
    })
    names(arl.gtypes) <- basename(arl.files)
    arl.gtypes
  }
  
  if(delete.files) file.remove(c(dir(label, full.names = T), label, file))
  
  result
}
