#' @rdname jack.hwe 
#' @usage \method{plot}{jack.influential}(x, main = "", ...)

plot.jack.influential <- function(x, main = "", ...) {
  odds.ratio <- as.vector(x$odds.ratio)
  odds.ratio <- odds.ratio[!is.na(odds.ratio) | !is.nan(odds.ratio) | !is.infinite(odds.ratio)]
  odds.sort <- sort(odds.ratio)
  odds.freq <- 1:length(odds.sort) / length(odds.sort)
  
  plot(odds.sort, odds.freq, type = "l", bty = "l", main = main,
       xlab = "Odds-ratio", ylab = "Cumulative frequency", ...
  )
  if(!all(is.null(x$influential))) {
    if(nrow(x$influential) > 0) {
      min.influential <- min(x$influential$odds.ratio)
      min.freq <- length(odds.ratio[odds.ratio < min.influential]) / length(odds.ratio)
      abline(v = min.influential, lty = "dashed", lwd = 1)    
      min.txt <- paste(round(min.influential, 2), " = ", round(100 * min.freq, 0), "%", sep = "")
      text(min.influential, min.freq, min.txt, adj = c(1.1, 1.4), cex = 1)        
    }
  }
}