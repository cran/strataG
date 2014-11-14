#' @rdname ti.tv.ratio

is.ti.tv <- function(b1, b2, sub.type) {
  if(!tolower(sub.type) %in% c("ti", "tv")) stop("'sub.type' must be 'ti' or 'tv'")
  b1 <- tolower(b1)
  b2 <- tolower(b2)
  if(!(all(c(b1, b2) %in% colnames(ti.tv.mat)))) return(FALSE)
  if(b1 == b2) return(FALSE)
  ti.tv.mat[b1, b2] == tolower(sub.type)
}