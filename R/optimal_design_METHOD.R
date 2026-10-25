# ========================== MÉTODOS DE SpatFDdesign ===========================

#' Print and plot a SpatFDdesign object
#'
#' `print()` shows the new stations and how the total prediction variance
#' decreased. `plot()` draws the candidates in grey, the existing stations in
#' black, the new stations in red, numbered in the order they were added, and,
#' when the targets were given as a map, the map.
#'
#' @param x An object of class `SpatFDdesign`, as returned by
#'   [optimal_design()].
#' @param main Title of the plot.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @name SpatFDdesign-print
#' @export
print.SpatFDdesign <- function(x, ...){
  cat('\n------------- Optimal design -------------\n')
  cat('Variable: ', attr(x, 'variable'), ' - Method: ', attr(x, 'method'), ' (ordinary kriging)\n', sep = '')
  cat('\nNew stations:\n')
  print(x$new_stations)
  cat('\nTotal prediction variance over the', nrow(attr(x, 'targets')), 'target points:\n')
  print(x$variance, row.names = FALSE, digits = 6)
  cat('------------------------------------------\n\n')
  invisible(x)
}

#' @rdname SpatFDdesign-print
#' @export
plot.SpatFDdesign <- function(x, main = NULL, ...){
  if (is.null(main)) main <- paste0('Optimal design: ', nrow(x$new_stations), ' new stations')

  stations   <- attr(x, 'stations')
  candidates <- attr(x, 'candidates')
  targets    <- attr(x, 'targets')
  map        <- attr(x, 'map')

  old <- par(ask = FALSE)
  on.exit(par(old), add = TRUE)

  everything <- rbind(stations, candidates, targets)
  plot(everything, type = 'n', asp = 1, xlab = colnames(stations)[1], ylab = colnames(stations)[2], main = main)

  # Sin mapa, los puntos objetivo se marcan; con mapa, el mapa los representa
  if (is.null(map)) {
    points(targets, pch = 3, cex = 0.5, col = 'grey70')
  } else {
    plot(sf::st_geometry(map), add = TRUE, border = 'grey50')
  }

  points(candidates, pch = 20, cex = 0.5, col = 'grey60')
  points(stations, pch = 19, col = 'black')
  points(x$new_stations, pch = 19, col = 'red', cex = 1.3)
  text(x$new_stations, labels = seq_len(nrow(x$new_stations)), pos = 3, col = 'red', font = 2)

  legend('topright', bty = 'n', cex = 0.85, pch = c(19, 19, 20), col = c('black', 'red', 'grey60'),
         legend = c('existing stations', 'new stations', 'candidates'))

  invisible(x)
}
