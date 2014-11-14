#' @export gelato.plot
#' 
#' @title Plot GELATo
#' @description Plot results from a GELATo run.
#' 
#' @param gelato.result the result of a call to \code{\link{gelato.run}}.
#' @param unknown the name of an unknown stratum in the \code{x$likelihoods} element.
#' @param main main label for top of plot.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{gelato.run}}

gelato.plot <- function(gelato.result, unknown, main = NULL) { 
  lik <- gelato.result$likelihoods[[unknown]]
  lik <- lik[!sapply(lik, is.null)]
  if(length(lik) == 0) stop(paste("No likelihood distributions available for '", unknown, "'", sep = ""))
  xticks <- pretty(unlist(sapply(lik, function(x) x$fst.dist)))
  xlim <- range(xticks)
  op <- par(mar = c(3, 3, 3, 2) + 0.1, oma = c(2, 2, 0.1, 0.1), mfrow = c(length(lik), 1))
  high.prob <- gelato.result$assign.prob[unknown, "assignment"]
  for(known in names(lik)) {
    known.lik <- lik[[known]]
    null.max <- max(hist(known.lik$fst.dist[, "null"], plot = F)$density)
    obs.max <- max(hist(known.lik$fst.dist[, "obs"], plot = F)$density)
    lik.mean <- known.lik$norm.coefs["mean"]
    lik.sd <- known.lik$norm.coefs["sd"]
    ylim <- range(pretty(c(0, null.max, obs.max, dnorm(lik.mean, lik.mean, lik.sd))))
    hist(known.lik$fst.dist[, "null"], breaks = 10, freq = FALSE, xlim = xlim, ylim = ylim, 
         xlab = "", ylab = "", main = "", col = "red", xaxt = "n")
    x <- NULL # To avoid R CMD CHECK warning about no global binding for 'x'
    curve(dnorm(x, lik.mean, lik.sd), from = xlim[1], to = xlim[2],
          add = TRUE, col = "black", lwd = 3, ylim = ylim)
    par(new = TRUE)
    hist(known.lik$fst.dist[, "obs"], breaks = 10, freq = FALSE, xlim = xlim, ylim = ylim, 
         xlab = "", ylab = "", col = "darkgreen", main = "", xaxt = "n", yaxt = "n") 
    axis(1, pretty(xlim))
    ll.median <- known.lik$log.Lik.smry["median"]
    log.lik <- if(!is.infinite(ll.median)) format(ll.median, digits = 4) else "Inf"
    p.val <- format(gelato.result$assign.prob[unknown, known], digits = 2)
    pop <- paste(known, " (lnL = ", log.lik, ", p = ", p.val, ")", sep = "")
    mtext(pop, side = 3, line = 1, adj = 1, font = ifelse(known == high.prob, 2, 1))
  }
  mtext("Fst", side = 1, outer = T, cex = 1.2)
  mtext("Density", side = 2, outer = T, cex = 1.2)
  par(op)
  if(!is.null(main)) mtext(main, side = 3, line = 3, adj = 0, font = 3)
}