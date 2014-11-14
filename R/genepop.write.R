#' @rdname genepop.run

genepop.write <- function(g, label = "genepop.write", input.fname = "loc_data.txt") {
  stopifnot.gtypes(g, "diploid")
  
  g <- decode(g)
  loc_dat <- rbind(sapply(1:length(attr(g, "locus.names")), function(i) {
    locus <- g$genotypes[, locus.cols(i), drop = FALSE]
    loc.levels <- sort(unique(as.vector(locus)))
    a1 <- as.numeric(factor(locus[, 1], levels = loc.levels))
    a2 <- as.numeric(factor(locus[, 2], levels = loc.levels))
    a1[is.na(a1)] <- 0
    a2[is.na(a2)] <- 0
    max.width <- max(2, nchar(a1), nchar(a2))
    a1 <- formatC(a1, width = max.width, flag = "0")
    a2 <- formatC(a2, width = max.width, flag = "0")
    paste(a1, a2, sep = "")
  }))
  
  id.vec <- sapply(rownames(g$genotypes), function(id) gsub(",", "_", id))
  pop.vec <- sapply(g$genotypes[, "strata"], function(pop) gsub(",", "_", pop))
  
  locus.names <- attr(g, "locus.names")
  names(locus.names) <- paste("LOC", 1:length(locus.names), sep = "")
  write(c(label, names(locus.names)), file = input.fname)
  for(pop in unique(pop.vec)) {
    write("POP", file = input.fname, append = TRUE)
    pop.rows <- which(pop.vec == pop)
    for(i in pop.rows) {
      write(paste(id.vec[i], pop, ",", paste(loc_dat[i, ], collapse = " ")), file = input.fname, append = TRUE)
    }  
  }
  
  invisible(locus.names)
}