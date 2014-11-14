#' @rdname read.arlequin

write.arlequin <- function(g, file = "gtypes.prj", title = "gtypes from R", data.type = "DNA") {
  stopifnot.gtypes(g)
  
  if(!data.type %in% c("DNA", "RFLP", "MICROSAT")) stop("'data.type' must be DNA, RFLP, or MICROSAT.")
  if(is.diploid(g) & data.type == "DNA") stop("'data.type' cannot be DNA if 'g' is diploid.")
  if(is.haploid(g) & data.type %in% c("RFLP", "MICROSAT")) stop("'data.type' cannot be RFLP or MICROSAT if 'g' is haploid.")
  
  write("[Profile]", file = file)
  write(paste("Title = \"", title, "\"", sep = ""), file = file, append = TRUE)
  write(paste("NbSamples =", nrow(g$genotypes)), file = file, append = TRUE)
  write(paste("DataType =", data.type), file = file, append = TRUE)
  write(paste("GenotypicData =", ifelse(is.haploid(g), 0, 1)), file = file, append = TRUE)
  
  g <- decode(g)
  write("[Data]", file = file, append = TRUE)
  if(!is.null(g$sequences)) {
    write("[[HaplotypeDefinition]]", file = file, append = TRUE)
    write("HaplListName=\"Haplotypes\"", file = file, append = TRUE)
    write("HaplList={", file = file, append = TRUE)
    for(x in names(g$sequences)) write(paste(x, g$sequences[[x]]), file = file, append = TRUE)
    write("}", file = file, append = TRUE)
  }
  
  invisible(NULL)
}