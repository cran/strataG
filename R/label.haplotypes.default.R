#' @rdname label.haplotypes
#' @importFrom swfscMisc zero.pad

label.haplotypes.default  <- function(x, prefix = NULL, ignore.gaps = FALSE) {  
  stopifnot.aligned(x)
  opt <- options(stringsAsFactors = FALSE)
  
  # make sure sequences have names
  if(is.null(names(x))) names(x) <- 1:length(x)
  
  # create matrix of aligned sequences and remove sites with n's
  seq.mat <- tolower(do.call(rbind, x))
  seq.mat <- seq.mat[, which(apply(seq.mat, 2, function(x) all(x != "n"))), drop = FALSE]
  
  hap.list <- if(all(seq.mat != "n") | nrow(seq.mat) == 1) {
    # if there are no N's do haplotype assignment based on simple grouping of unique sequences (faster)
    if(ignore.gaps) {
      no.gaps <- apply(seq.mat, 2, function(x) all(x != "-"))
      seq.mat <- seq.mat[, no.gaps, drop = FALSE]
    }
    seq.vec <- apply(seq.mat, 1, paste, collapse = "")
    result <- lapply(unique(seq.vec), function(x) names(seq.vec)[seq.vec == x])
    names(result) <- sapply(result, function(h) h[1])
    result
  } else {
    # otherwise do pairwise comparison of sequences and identify which pairs are the same
    pair.same.mat <- matrix(as.logical(NA), nrow = nrow(seq.mat), ncol = nrow(seq.mat))
    diag(pair.same.mat) <- TRUE
    rownames(pair.same.mat) <- colnames(pair.same.mat) <- rownames(seq.mat)
    pairs <- combn(rownames(seq.mat), 2) 
    for(i in 1:ncol(pairs)) {
      id1 <- pairs[1, i]
      id2 <- pairs[2, i]
      pair.seq <- seq.mat[c(id1, id2), , drop = FALSE]
      to.keep <- which(apply(pair.seq, 2, function(bases) {
        no.gaps <- if(ignore.gaps) "-" %in% bases else TRUE
        no.gaps & all(bases != "n")
      }))
      pair.seq <- pair.seq[, to.keep, drop = FALSE]
      pair.same.mat[id1, id2] <- pair.same.mat[id2, id1] <- all(pair.seq[1, ] == pair.seq[2, ])
    }
    
    # find all ids shared by each sample
    result <- if(sum(pair.same.mat) > nrow(pair.same.mat)) {
      lapply(rownames(pair.same.mat), function(i) {
        id.row <- pair.same.mat[i, , drop = TRUE]
        sort(names(id.row)[id.row])
      })
    } else as.list(rownames(seq.mat))
    # rename haplotype list with first id found in shared list
    names(result) <- sapply(result, function(h) h[1])
    # remove duplicate names
    result[!duplicated(names(result))]
  }
  
  # sort by frequency
  hap.list <- hap.list[order(sapply(hap.list, length), decreasing = TRUE)]
  
  # rename haplotype list if prefix given
  if(!is.null(prefix)) {
    hap.nums <- zero.pad(1:length(hap.list))  
    names(hap.list) <- paste(prefix, hap.nums, sep = "")
  }
  
  # get haplotype sequences
  hap.seqs <- lapply(hap.list, function(i) tolower(create.consensus(x[i], ignore.gaps = ignore.gaps)))
  names(hap.seqs) <- names(hap.list)
  
  # create haplotype assignment vector
  hap.df <- do.call(rbind, lapply(names(hap.list), function(i) data.frame(id = hap.list[[i]], haplotype = i)))
  hap.df <- hap.df[match(names(x), hap.df$id), ]
  haps <- hap.df$haplotype
  names(haps) <- hap.df$id
  
  options(opt)
  list(haps = haps, hap.seqs = hap.seqs)
}