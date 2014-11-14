#' @rdname phase.run

phase.posterior <- function(ph.res, keep.missing = TRUE) {    
  if(!"phase.result" %in% class(ph.res)) stop("'ph.res' is not a result from 'phase.run'.")
                                                                    
  num.iter <- length(ph.res[[1]]$posterior)
  lapply(1:num.iter, function(iter) {
    do.call(cbind, lapply(1:length(ph.res), function(locus) {
      ph.res <- ph.res[[locus]]
      post.df <- ph.res$posterior[[iter]]
      
      if(keep.missing) {
        for(i in 1:nrow(post.df)) {
          row <- which(rownames(ph.res$orig.gtypes$genotypes) == post.df[i, 1])
          if(any(ph.res$orig.gtypes$genotypes[row, -1] == -1)) post.df[i, 2:3] <- NA
        }
      }
      
      cols <- if(locus == 1) {1:3} else {2:3}
      post.df[, cols]
    })) 
  })
}