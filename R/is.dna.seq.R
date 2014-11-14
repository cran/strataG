#' @rdname as.dna.seq

is.dna.seq <- function(x, show.warnings = FALSE) {  
  if(is.null(x)) return(FALSE)
  if(is.null(names(x)) & show.warnings) warning("no sample ids are associated")
  if(!is.list(x)) {
    if(show.warnings) warning("'x' should be a list")
    return(FALSE)
  }
  if(!all(sapply(x, function(i) is.vector(i) & is.character(i)))) {
    if(show.warnings) warning("each element of 'x' should be a character vector")
    return(FALSE)
  }
  one.char <- all(sapply(x, function(i) all(sapply(i, nchar)) == 1)) 
  if(!one.char) {
    if(show.warnings) warning("each element of each sequence in 'x' must be one character")
    return(FALSE)
  }
  invalid.chars <- any(sapply(x, function(i) !all(tolower(i) %in% rownames(iupac.mat)))) 
  if(invalid.chars) {
    if(show.warnings) warning("'x' includes invalid characters")
    return(FALSE)
  }
  TRUE
}