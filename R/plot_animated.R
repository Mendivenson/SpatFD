# =============================== plot_animated() ==============================

#' Interactive maps of a prediction through time
#'
#' Draws, with \pkg{plotly}, the maps of a prediction made over a map: on the
#' left the predicted curves at one time, with a slider (and a play button) to
#' move through time; on the right the prediction variance, which does not
#' depend on time. In RStudio it opens in the Viewer; in an R Markdown HTML
#' document it stays interactive.
#'
#' @param x An object of class `SpatFDpred`, as returned by
#'   [predict()][predict.SpatFDmodel] with a map as `newcoords`.
#' @param variable Name of the variable to plot. Defaults to the first
#'   predicted variable.
#' @param times Times of the slider. Either a single number, the number of
#'   equally spaced times between the first and the last time of the curves,
#'   or a vector of times in the same format as the times of the data
#'   (numbers, dates or dates and times).
#' @param main Title of the plot. Defaults to the name of the variable.
#'
#' @details
#' The colour scale of the prediction is the same at every time, so the
#' colours can be compared between times. The observed sites are drawn as
#' white points. Every time of the slider adds a copy of the prediction map
#' to the plot, so many times or a fine grid make a large plot.
#'
#' Requires the packages \pkg{plotly} and \pkg{sf}.
#'
#' @return A \pkg{plotly} object (an HTML widget).
#'
#' @seealso [plot()][plot.SpatFDpred] for static maps.
#'
#' @examplesIf requireNamespace("plotly", quietly = TRUE) && requireNamespace("sf", quietly = TRUE)
#' set.seed(1)
#' dat <- sapply(1:10, function(i) rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1))
#' colnames(dat) <- paste0("s", 1:10)
#' coords <- cbind(x = runif(10), y = runif(10))
#' rownames(coords) <- colnames(dat)
#' sim <- spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "sim")
#'
#' table <- data.frame(var = c("sim.PC1", "sim.PC2"), model = "exponential",
#'                     rango = 0.5, sigma = c(30, 0.1))
#' model  <- spatfd_model(sim, model = table)
#' square <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))))
#'
#' plot_animated(predict(model, square, grid_size = 400), times = 20)
#'
#' @export
plot_animated <- function(x,                 # Predicción sobre un mapa, ver predict()
                          variable = NULL,   # Variable que se grafica
                          times = 100,       # Número de tiempos o tiempos del slider
                          main = NULL){      # Título

  if (!requireNamespace('plotly', quietly = TRUE) || !requireNamespace('sf', quietly = TRUE)) {
    stop("plot_animated() needs the packages 'plotly' and 'sf'", call. = F)
  }
  stopifnot("'x' must be a prediction, see predict()" = inherits(x, 'SpatFDpred'))
  if (is.null(attr(x, 'map'))) {
    stop('plot_animated() is only available for predictions made over a map', call. = F)
  }

  if (is.null(variable)) variable <- names(x)[1]
  if (!(variable %in% names(x))) {
    stop('Variable ', variable, ' was not predicted. The predicted variables are: ',
         paste(names(x), collapse = ', '), call. = F)
  }
  if (is.null(main)) main <- variable

  internal <- attr(x, 'internal')
  curves   <- internal$pred[[variable]]
  days     <- slider_days(times, curves$basis$rangeval, internal$time[[variable]])
  labels   <- slider_labels(days, internal$time[[variable]])

  # VALORES: La predicción en cada tiempo del slider y la varianza, como
  # matrices sobre la grilla. Se guardan con 4 cifras significativas, que no
  # cambian los colores y hacen el gráfico mucho más liviano.
  cells      <- grid_cells(attr(x, 'newcoords'))
  predicted  <- signif(fda::eval.fd(days, curves), 4)
  prediction <- lapply(seq_along(days), function(k) cells$fill(predicted[k, ]))
  variance   <- cells$fill(signif(x[[variable]]$var, 4))

  animated_maps(cells, prediction, variance, labels, map_border(attr(x, 'map')),
                internal$coords[[variable]], range(predicted), main)
}
