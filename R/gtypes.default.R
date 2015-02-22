#' @rdname gtypes

gtypes.default <- function(gen.data, id.col = 1, strata.col = 2, locus.col = 3, dna.seq = NULL, description = NULL, delete.missing.strata = TRUE, code.start = 0, ...) {
  # Check that gen.data is a matrix or data.frame
  if(!(is.matrix(gen.data) | is.data.frame(gen.data))) stop("'gen.data' must be a matrix or data.frame")
  
  if(length(id.col) > 1 & !(is.character(id.col) | is.numeric(id.col))) stop("'id.col' must be a one-element character or numeric vector")
  if(length(strata.col) > 1 & !(is.character(strata.col) | is.numeric(strata.col))) stop("'strata.col' must be a one-element character or numeric vector")
  if(length(locus.col) > 1 & !(is.character(locus.col) | is.numeric(locus.col))) stop("'locus.col' must be a one-element numeric vector")
  
  if(is.character(id.col)) if(!id.col %in% colnames(gen.data)) stop(paste("'", id.col, "' is not a column in 'gen.data'", sep = ""))
  if(is.character(strata.col)) if(!strata.col %in% colnames(gen.data)) stop(paste("'", strata.col, "' is not a column in 'gen.data'", sep = ""))
  if(locus.col > ncol(gen.data)) stop("'locus.col' can't be greater than the number of columns in 'gen.data'")
    
  # Extract genetic data based on specified columns
  gen.data <- cbind(gen.data[, id.col, drop = FALSE], gen.data[, strata.col, drop = FALSE], gen.data[, locus.col:ncol(gen.data), drop = FALSE])
  gen.data.cols <- colnames(gen.data)
  gen.data <- do.call(cbind, lapply(1:ncol(gen.data), function(i) as.character(gen.data[, i])))
  colnames(gen.data) <- gen.data.cols
  gen.data <- gen.data[!is.na(gen.data[, 1]), , drop = FALSE]
  if(delete.missing.strata) gen.data <- gen.data[!is.na(gen.data[, 2]), , drop = FALSE]
  
  ids <- gen.data[, 1]
  id.freq <- table(ids)
  if(any(id.freq > 1)) {
    dup.ids <- paste(names(id.freq)[id.freq > 1], collapse = ", ")
    stop(paste("The following IDs in 'locus.data' occur more than once:", dup.ids))
  }
  
  strata <- gen.data[, 2]
  names(strata) <- ids
  
  locus.data <- gen.data[, 3:ncol(gen.data), drop = FALSE]
  rownames(locus.data) <- ids
  is.haploid <- ncol(locus.data) == 1
  
  # Extract locus.names
  locus.names <- locus.names(locus.data)
  colnames(locus.data) <- if(ncol(locus.data) > 1) { 
    paste(rep(locus.names, each = 2), 1:2, sep = ".")
  } else {
    locus.names
  } 

  # Check for loci with one or both alleles missing for all samples
  if(ncol(locus.data) > 1) {
    missing.data.loci <- which(sapply(1:num.loci(locus.data), function(i) {
      this.loc <- locus.data[, locus.cols(i) - 1]
      all(is.na(this.loc[, 1])) | all(is.na(this.loc[, 2]))
    }))
    if(length(missing.data.loci) > 0) {
      missing.names <- paste(locus.names[missing.data.loci], collapse = ", ")
      warning(paste("The following loci will be deleted because they are missing all data for one or both alleles:", missing.names))
      to.delete <- as.vector(sapply(missing.data.loci, function(i) locus.cols(i) - 1))
      locus.data <- locus.data[, -to.delete]
      locus.names <- locus.names(locus.data)
    }  
  }
  
  # Check dna.seq
  if(!is.null(dna.seq)) {
    if(!is.dna.seq(dna.seq)) stop("'dna.seq' is not a valid list of sequences")
    if(!is.haploid) stop("Can't include 'dna.seq' with diploid data")
    haps <- unique(locus.data[, 1])
    missing.haps <- haps[!(haps %in% names(dna.seq))]
    if(length(missing.haps) > 0) {
      missing.haps <- paste(missing.haps, collapse = ", ")
      stop(paste("The following haplotype labels were not found in 'dna.seq':", missing.haps))
    }
    dna.seq <- lapply(dna.seq, tolower)
  }
  
  # Encode strata, locus.data, and sequence names by replacing characters with numerics
  
  convert.fac <- function(x.fac, start) {
    nums <- as.numeric(x.fac) - 1 + start
    nums[is.na(nums)] <- -1
    num.names <- c(NA, levels(x.fac))
    names(num.names) <- c("-1", 1:(length(num.names) - 1) - 1 + start)
    list(nums = nums, num.names = num.names)
  }
  
  # Convert strata
  strata.nums <- convert.fac(factor(strata), code.start)
  recoded.strata <- strata.nums$nums
  names(recoded.strata) <- names(strata)
  strata.names <- strata.nums$num.names
  
  # Convert locus.data
  recoded.loci <- if(ncol(locus.data) == 1) {
    haps <- convert.fac(factor(locus.data[, 1]), code.start)
    result <- matrix(haps$nums, ncol = 1)
    list(loci = result, allele.names = list(haps$num.names))
  } else {
    locus.list <- lapply(seq(1, ncol(locus.data), 2), function(i) {
      locus <- locus.data[, c(i, i + 1)]
      alleles <- sort(unique(as.vector(locus)))
      convert.fac(factor(locus, levels = alleles), code.start)
    })
    list(loci = do.call(cbind, lapply(locus.list, function(x) matrix(x$nums, ncol = 2))),
         allele.names = lapply(locus.list, function(x) x$num.names)
    )
  }
  rownames(recoded.loci$loci) <- rownames(locus.data)
  colnames(recoded.loci$loci) <- colnames(locus.data)
  allele.names <- recoded.loci$allele.names
  names(allele.names) <- locus.names
  
  # Convert names of sequences
  if(!is.null(dna.seq)) {
    dna.seq <- dna.seq[allele.names[[1]][-1]]
    names(dna.seq) <- names(allele.names[[1]])[-1]
  }
  
  recoded.strata <- recoded.strata[rownames(recoded.loci$loci)]
  g <- list(genotypes = cbind(strata = recoded.strata, recoded.loci$loci), sequences = dna.seq)
  attr(g, "strata.names") <- strata.names
  attr(g, "locus.names") <- locus.names
  attr(g, "allele.names") <- recoded.loci$allele.names
  attr(g, "description") <- if(is.null(description)) paste("gtypes created on", format(Sys.time(), "%Y-%m-%d %H:%M")) else description
  class(g) <- c("gtypes", class(g))
  options(check.gtypes = TRUE)
    
  # Remove samples with missing data for all loci
  na.omit(g, strata = delete.missing.strata, loci = TRUE)
}