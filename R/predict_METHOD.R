# ========================== MÉTODOS DE SpatFDpred =============================
# Impresión y gráficos de las predicciones. Las curvas funcionales (predichas y
# observadas) viven en el atributo oculto 'internal' de la predicción.

#' Print a SpatFDpred object
#'
#' Shows the number of prediction sites, the predicted variables, the method
#' and the prediction variance of the first sites.
#'
#' @param x An object of class `SpatFDpred`, as returned by
#'   [predict()][predict.SpatFDmodel].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.SpatFDpred <- function(x, ...){
  sites <- nrow(attr(x, 'newcoords'))
  shown <- min(sites, 10)

  cat('\n\n-------------------------------------------------\n')
  cat('Predicted curves at ', sites, ' sites', if (!is.null(attr(x, 'map'))) ' over a map',
      ' for: ', paste(names(x), collapse = ', '), '\n', sep = '')
  cat('Method: ', attr(x, 'method'), ', ordinary kriging\n', sep = '')
  cat('\nPrediction variance:\n')
  print(do.call(cbind, lapply(x, function(variable) variable$var[seq_len(shown)])))
  if (shown < sites) cat('... and', sites - shown, 'more sites\n')
  cat('-------------------------------------------------\n\n')
  invisible(x)
}

#' Plot the predictions of a SpatFDmodel object
#'
#' Plots the curves predicted by [predict()][predict.SpatFDmodel] for one
#' variable, either on top of the observed curves or, when the prediction was
#' made over a map, as two maps side by side: the prediction and its variance.
#'
#' @param x An object of class `SpatFDpred`, as returned by
#'   [predict()][predict.SpatFDmodel].
#' @param variable Name of the variable to plot. Defaults to the first
#'   predicted variable.
#' @param mode Either `"curves"`, to draw the predicted curves over the
#'   observed ones, or `"map"`, to draw the maps of the prediction and of its
#'   variance. Defaults to `"map"` when the prediction was made over a map and
#'   to `"curves"` otherwise.
#' @param time Only for `mode = "map"`: time at which the predicted curves are
#'   shown in the prediction map, in the same format as the times of the data
#'   (a number, a date or a date and time). Defaults to the first time.
#' @param main Title of the plot.
#' @param ... Graphical parameters. For `mode = "curves"` they are applied to
#'   the predicted curves (for example `col` or `lwd`); for `mode = "map"`
#'   only `col`, the colour palette of the maps, is used.
#'
#' @details
#' With `mode = "curves"` the observed curves are drawn in grey and the
#' predicted ones in blue. When the times of the data are dates, the
#' horizontal axis shows them as dates.
#'
#' A predicted curve is a function of time, so a map can only show it at one
#' time: the one given in `time`. The variance map does not depend on `time`,
#' because the prediction variance is that of the whole curve. The observed
#' sites are drawn as white points.
#'
#' @return `x`, invisibly. Called for its side effect of drawing a plot.
#'
#' @seealso [predict()][predict.SpatFDmodel]
#'
#' @examples
#' set.seed(1)
#' dat <- sapply(1:10, function(i) rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1))
#' colnames(dat) <- paste0("s", 1:10)
#' coords <- cbind(x = runif(10), y = runif(10))
#' rownames(coords) <- colnames(dat)
#' sim <- spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "sim")
#'
#' table <- data.frame(var = c("sim.PC1", "sim.PC2"), model = "exponential",
#'                     rango = 0.5, sigma = c(30, 0.1))
#' model <- spatfd_model(sim, model = table)
#'
#' # Predicted curves over the observed ones
#' plot(predict(model, newcoords = cbind(x = c(0.2, 0.7), y = c(0.5, 0.5))))
#'
#' # Prediction over a map: maps of the prediction and of its variance
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   square <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))))
#'   plot(predict(model, newcoords = square, grid_size = 400), time = 25)
#' }
#'
#' @export
plot.SpatFDpred <- function(x, variable = NULL, mode = NULL, time = NULL, main = NULL, ...){

  if (is.null(variable)) {
    variable <- names(x)[1]
    if (length(x) > 1) message('Plotting the variable ', variable)
  }
  if (!(variable %in% names(x))) {
    stop('Variable ', variable, ' was not predicted. The predicted variables are: ',
         paste(names(x), collapse = ', '), call. = F)
  }

  if (is.null(mode)) mode <- if (is.null(attr(x, 'map'))) 'curves' else 'map'
  if (is.null(main)) main <- variable

  # El gráfico no pide <Enter>, aunque el dispositivo esté en modo 'ask'
  old <- par(ask = FALSE)
  on.exit(par(old), add = TRUE)

  switch(mode,
         'curves' = plot_curves(x, variable, main, ...),
         'map'    = plot_maps(x, variable, time, main, ...),
         stop("'mode' must be the string 'curves' or the string 'map'", call. = F))

  invisible(x)
}

# CURVAS PREDICHAS SOBRE LAS OBSERVADAS: En gris las curvas observadas
# (suavizadas) y en azul las predichas, con el eje del tiempo en el formato de
# los datos (fechas, si eran fechas)
plot_curves <- function(x, variable, main, ...){
  internal  <- attr(x, 'internal')
  observed  <- internal$observed[[variable]]
  predicted <- internal$pred[[variable]]
  template  <- internal$time[[variable]]

  grid <- function(fd) seq(fd$basis$rangeval[1], fd$basis$rangeval[2], length.out = 400)
  t_obs  <- grid(observed)
  t_pred <- grid(predicted)
  y_obs  <- fda::eval.fd(t_obs, observed)
  y_pred <- fda::eval.fd(t_pred, predicted)

  style <- utils::modifyList(list(col = 'blue', lty = 1, lwd = 1.5), list(...))

  # El eje se arma con los tiempos del usuario; las líneas usan su valor
  # numérico interno (días para Date, segundos para POSIXct)
  axis_obs  <- from_days(t_obs, template)
  axis_pred <- from_days(t_pred, template)

  plot(range(axis_obs), range(y_obs, y_pred), type = 'n', xlab = 'time', ylab = variable, main = main)
  matlines(as.numeric(axis_obs), y_obs, col = 'grey80', lty = 1)
  do.call(matlines, c(list(as.numeric(axis_pred), y_pred), style))
}

# MAPAS DE PREDICCIÓN Y VARIANZA, LADO A LADO: La predicción se muestra en un
# tiempo dado, la varianza es la de toda la curva. De arriba hacia abajo:
#
#     +------------------- título -------------------+
#     +------------ subtítulo con detalles ----------+
#     |  barra de color     |    barra de color      |
#     |  mapa               |    mapa                |
#     |  Prediction         |    Variance            |
#     +----------------------------------------------+
plot_maps <- function(x, variable, time, main, ...){
  map <- attr(x, 'map')
  if (is.null(map)) stop("mode = 'map' is only available for predictions made over a map", call. = F)

  internal <- attr(x, 'internal')
  curves   <- internal$pred[[variable]]
  times    <- curves$basis$rangeval

  # El tiempo se da en el formato de los datos y se pasa a días, como los
  # tiempos de las curvas
  if (is.null(time)) time <- from_days(times[1], internal$time[[variable]])
  day <- as_days(time)
  stopifnot(
    "'time' must be a single value within the times of the curves" =
      is.numeric(day) && length(day) == 1 && day >= times[1] && day <= times[2]
  )

  details  <- paste0('Time: ', format(time), ' - Method: ', attr(x, 'method'), ' (ordinary kriging)')
  grid     <- attr(x, 'newcoords')
  stations <- internal$coords[[variable]]
  colours  <- utils::modifyList(list(col = hcl.colors(100, 'viridis', rev = TRUE)), list(...))$col

  # Los puntos forman una grilla regular: se pasan a una matriz para image()
  cells <- grid_cells(grid)

  # Cada mapa ocupa dos paneles: su barra de color y, debajo, el mapa
  draw <- function(values, what){
    limits <- range(values)
    if (diff(limits) == 0) limits <- limits + c(-0.5, 0.5)

    par(mar = c(2, 6, 0.5, 6))
    image(x = seq(limits[1], limits[2], length.out = length(colours)), y = 1,
          z = matrix(seq_along(colours), ncol = 1), col = colours,
          axes = FALSE, xlab = '', ylab = '')
    axis(1, cex.axis = 1.3, mgp = c(3, 0.6, 0), tcl = -0.3)
    box()

    par(mar = c(3, 1, 0.5, 1))
    image(x = cells$x, y = cells$y, z = t(cells$fill(values)), zlim = limits, col = colours, asp = 1,
          axes = FALSE, xlab = '', ylab = '')
    plot(sf::st_geometry(map), add = TRUE)
    points(stations, pch = 21, bg = 'white', cex = 1.2)

    mtext(what, side = 1, line = 1.2, cex = 1.3, font = 2)
  }

  old <- par(mar = c(0, 0, 0, 0))
  on.exit({ layout(1); par(old) })

  layout(matrix(c(1, 1,
                  2, 4,
                  3, 5), nrow = 3, byrow = TRUE), heights = c(2, 1.3, 10))

  plot.new()
  text(0.5, 0.68, main, cex = 2.4, font = 2)
  text(0.5, 0.18, details, cex = 1.5)

  draw(as.vector(fda::eval.fd(day, curves)), 'Prediction')
  draw(x[[variable]]$var, 'Variance')
}
