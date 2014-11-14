#' @rdname read.mega

write.mega <- function(g, file, title = NULL, line.width = 60) {
  stopifnot.gtypes(g)
  
  if(is.null(title)) title <- attr(g, "description")
  dna.seq <- decode.sequences(g)
  
  write("#MEGA", file)
  write(paste("title:", title, sep = ""), file, append = TRUE)
  write("", file, append = TRUE)
  for(x in names(dna.seq)) {
    write(paste("#", x, sep = ""), file, append = TRUE)
    mt.seq <- dna.seq[[x]]
    for(j in seq(1, length(mt.seq), by = line.width)) {
      seq.line <- paste(mt.seq[j:(j + line.width - 1)], collapse = "")
      write(seq.line, file, append = TRUE)
    }
    write("", file, append = TRUE)
  }
}