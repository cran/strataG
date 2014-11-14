#' @rdname phase.run

phase.write <- function(g, loci, positions = NULL, type = rep("S", length(loci)), in.file = "phase_in") {  
  stopifnot.gtypes(g, "diploid")
  
  # Make sure locus.names and locus.positions are sorted properly
  if(is.null(positions)) positions <- rep(1, length(attr(g, "locus.names")))
  asc.order <- order(positions)
  loci <- loci[asc.order]
  positions <- positions[asc.order]

  sub.g <- subset(g, loci = loci)
  g <- decode(sub.g)
  write(c(
    nrow(g$genotypes), 
    length(loci),
    paste("P", paste(positions, collapse = " ")),
    paste(type, collapse = ""), ""
  ), file = in.file)
  
  g$genotypes[is.na(g$genotypes)] <- "?"
  for(i in 1:nrow(g$genotypes)) {
    write(c(
      rownames(g$genotypes)[i],
      paste(g$genotypes[i, seq(2, ncol(g$genotypes) - 1, 2)], collapse = " "),
      paste(g$genotypes[i, seq(3, ncol(g$genotypes), 2)], collapse = " ")
    ), file = in.file, append = TRUE)
  }
  
  invisible(list(filename = in.file, gtypes = sub.g))
}