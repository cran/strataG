#' @rdname locus.names

locus.cols <- function(x, g = NULL) {
 if(is.character(x) | is.numeric(x)) {
    if(is.character(x)) {
      if(!is.gtypes(g)) stop("if 'x' is character, then 'g' must be a valid gtypes object.")
      if(!all(x %in% attr(g, "locus.names"))) stop("some 'x' cannot be found in 'g'.")
    } else if(!all(x > 0)) stop("if numeric, 'x' must be positive integer(s)")
  } else stop("'x' must be character or numeric vector")

  sapply(x, function(i) {
    if(is.character(i)) i <- match(i, attr(g, "locus.names"))
    c((i * 2) - 1, i * 2) + 1
  })
}