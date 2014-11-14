#' @rdname read.mega
#' @export read.mega write.mega
#' @aliases mega
#' 
#' @title Read and Write MEGA
#' @description Read and write MEGA formatted files.
#' 
#' @param file a MEGA-formatted file of sequences.
#' @param g a \code{\link{gtypes}} object.
#' @param title title for data in file.
#' @param line.width width of sequence lines.
#' 
#' @return for \code{read.mega}, a list of:
#' \tabular{ll}{
#'   \code{title} \tab title of MEGA file.\cr
#'   \code{dna.seq} \tab a list of DNA sequences.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}

read.mega <- function(file) {  
  mega.file <- scan(file, what = "character", sep = "\n", strip.white = TRUE, quiet = TRUE)
  markers <- c(grep("#", mega.file), length(mega.file) + 1)
  title <- paste(mega.file[2:(markers[2] - 1)], collapse = " ")
  seq.df <- as.data.frame(t(sapply(2:(length(markers) - 1), function(i) {
    id <- sub("#", "", mega.file[markers[i]])
    seq.start <- markers[i] + 1
    seq.stop <- markers[i + 1] - 1
    dna.seq <- paste(mega.file[seq.start:seq.stop], collapse = "")
    c(id = id, sequence = dna.seq)
  })), stringsAsFactors = FALSE)
  dna.seq <- strsplit(seq.df$dna.seq, "")
  names(dna.seq) <- dna.seq$id
  list(title = title, dna.seq = dna.seq)
}