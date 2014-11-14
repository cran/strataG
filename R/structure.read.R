#' @rdname structure.run

structure.read <- function(file, pops = NULL) {
  if(!file.exists(file)) stop("'file' can't be found.")
  
  # Read file to get single results and parameter values
  result <- scan(file, "character", quiet = TRUE)
  
  loc <- grep("Estimated", result, ignore.case = FALSE, value = FALSE)
  est.ln.prob <- as.numeric(result[loc[1] + 6])
  
  loc <- grep("likelihood", result, ignore.case = FALSE, value = FALSE)
  mean.lnL <- as.numeric(result[loc[1] + 2])
  var.lnL <- as.numeric(result[loc[2] + 2])
   
  loc <- grep("MAXPOPS", result, value = F)
  maxpops <- result[loc]
  maxpops <- sub("MAXPOPS=", "", maxpops)
  maxpops <- as.integer(sub(",", "", maxpops))
  
  loc <- grep("GENSBACK", result, value = F)
  gensback <- result[loc]
  gensback <- sub("GENSBACK=", "", gensback)
  gensback <- as.integer(sub(",", "", gensback))
  
  smry <- c(k = maxpops, est.ln.prob = est.ln.prob, mean.lnL = mean.lnL, var.lnL = var.lnL)
  
  # Read file to get population assignment probability table
  result <- scan(file, "character", sep = "\n", quiet = TRUE)
  first <- grep("(%Miss)", result, value = FALSE) + 1
  last <- grep("Estimated Allele", result, value = FALSE) - 1
  tbl.txt <- result[first:last]
  
  # Remove special characters from table and whitespace from end of lines
  tbl.txt <- sub("[*]+", "", tbl.txt)
  tbl.txt <- sub("[(]", "", tbl.txt)
  tbl.txt <- sub("[)]", "", tbl.txt)
  tbl.txt <- sub("[|][ ]+$", "", tbl.txt)
  
  # Find which lines have had population priors
  prior.lines <- grep("[|]", tbl.txt)
  
  # Create table of population assignments for lines without population priors
  no.prior <- if(length(prior.lines) < length(tbl.txt)) {
    no.prior.q.txt <- if(length(prior.lines) == 0) tbl.txt else tbl.txt[-prior.lines]
    structure.parse.q.mat(no.prior.q.txt, pops)
  } else NULL
  
  # Return just this base table if MAXPOPS = 1
  if(maxpops == 1) {
    no.prior$row <- NULL
    return(list(summary = smry, q.mat = no.prior, prior.anc = NULL))
  }
  
  # Create table of population assignments for lines with population priors
  has.prior <- if(length(prior.lines) > 0) {
    prior.txt <- strsplit(tbl.txt[prior.lines], "[|]")
    # Get base table
    prior.q.txt <- unlist(lapply(prior.txt, function(x) x[1]))
    df <- structure.parse.q.mat(prior.q.txt, pops)
    # Parse ancestry assignments into matrix
    prior.anc <- lapply(prior.txt, function(x) {
      anc.mat <- matrix(NA, nrow = maxpops, ncol = gensback + 1)
      rownames(anc.mat) <- paste("Pop", 1:nrow(anc.mat), sep = ".")
      colnames(anc.mat) <- paste("Gen", 0:gensback, sep = ".")
      # Split on whitespace and colons
      x <- sapply(strsplit(x[-1], "\\s|[:]"), function(y) {
        y <- y[y != ""] # remove empty strings
        y[-1] # return vector with first element ("Pop") removed - vector has population # and ancestry assignments
      })
      # Populate ancestry matrix with probabilities
      for(i in 1:ncol(x)) {
        pop <- as.numeric(x[1, i])
        anc.mat[pop, ] <- as.numeric(x[-1, i])
      }
      anc.mat
    })
    names(prior.anc) <- df$id
    
    # Create population probability matrix for samples with priors and add to end of base table
    prob.mat <- t(sapply(1:nrow(df), function(i) {
      pop.probs <- rowSums(prior.anc[[i]])
      pop.probs[is.na(pop.probs)] <- df$prob.1[i]
      pop.probs
    }))
    colnames(prob.mat) <- paste("prob", 1:ncol(prob.mat), sep = ".")
    df$prob.1 <- NULL
    df <- cbind(df, prob.mat)
    
    list(df = df, prior.anc = prior.anc)
  } else NULL
  
  # Combine assignment probability matrices
  has.prior.df <- if(is.null(has.prior)) NULL else has.prior$df
  q.mat <- rbind(no.prior, has.prior.df)
  q.mat <- q.mat[order(q.mat$row), ]
  q.mat$row <- NULL
  rownames(q.mat) <- NULL
  # Make sure all probs sum to 1
  q.mat[, -(1:3)] <- t(apply(q.mat[, -(1:3)], 1, function(i) i / sum(i)))
  
  prior.anc <- if(is.null(has.prior)) NULL else has.prior$prior.anc
  
  list(summary = smry, q.mat = q.mat, prior.anc = prior.anc)
}