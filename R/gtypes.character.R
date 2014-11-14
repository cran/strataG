#' @rdname gtypes

gtypes.character <- function(gen.data, ...) {
  gen.data.mat <- if(length(gen.data) == 1) {
    if(!file.exists(gen.data)) stop(paste("A file named, '", gen.data, "' cannot be found", sep = ""))
    read.gen.data(gen.data)
  } else gen.data
  gtypes.default(gen.data.mat, ...)  
}
