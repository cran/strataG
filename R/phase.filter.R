#' @rdname phase.run

phase.filter <- function(ph.res, thresh = 0.5, keep.missing = TRUE) {
  if(!"phase.result" %in% class(ph.res)) stop("'ph.res' is not a result from 'phase.run'.")
  
  filtered <- lapply(ph.res, function(x) {
    gtype.probs <- x$gtype.probs
    locus.filtered <- do.call(rbind, lapply(unique(gtype.probs[, 1]), function(i) {
      this.id <- gtype.probs[gtype.probs[, 1] == i, ]
      max.index <- which.max(this.id$pr)
      if(length(max.index) == 0) return(this.id[1, ])
      kept.line <- this.id[max.index, ]
      if(as.numeric(kept.line$pr) < thresh) kept.line[, 2:3] <- c(NA, NA)
      kept.line
    })) 
    rownames(locus.filtered) <- NULL
    
    if(keep.missing) {
      for(i in 1:nrow(locus.filtered)) {
        row <- which(rownames(x$orig.gtypes$genotypes) == locus.filtered[i, 1])
        if(any(x$orig.gtypes$genotypes[row, -1] == -1)) locus.filtered[i, 2:3] <- NA
      }
    }
    
    locus.filtered
  })
  
  id <- filtered[[1]][, 1]
  filtered <- as.matrix(do.call(cbind, lapply(filtered, function(x) x[, 2:3])))
  rownames(filtered) <- id
  colnames(filtered) <- paste(rep(names(ph.res), each = 2), ".", 1:2, sep = "")
  filtered
}