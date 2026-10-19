# ========================== MÉTODOS DE SpatFDpred =============================
# Impresión y gráficos de las predicciones. Las curvas funcionales (predichas y
# observadas) viven en el atributo oculto 'internal' de la predicción.

#' @rdname predict.SpatFD
#' @param x An object of class `SpatFDpred`, as returned by `predict()`.
#' @export
print.SpatFDpred <- function(x, ...){
  sites <- nrow(attr(x, 'newcoords'))
  shown <- min(sites, 10)

  cat('\n\n-------------------------------------------------\n')
  cat('Predicted curves at ', sites, ' sites', if (!is.null(attr(x, 'map'))) ' over a map',
      ' for: ', paste(names(x), collapse = ', '), '\n', sep = '')
  cat('Method: ', attr(x, 'method'), ', ', attr(x, 'kriging'), ' kriging',
      if (isTRUE(attr(x, 'crossed'))) ', with cross-variograms', '\n', sep = '')
  cat('\nPrediction variance:\n')
  print(do.call(cbind, lapply(x, function(variable) variable$var[seq_len(shown)])))
  if (shown < sites) cat('... and', sites - shown, 'more sites\n')
  cat('-------------------------------------------------\n\n')
  invisible(x)
}

# FUNCIÓN DE PLOTEO: La función plot.SpatFDpred() tiene dos modos:
#                      'curves' -> Curvas predichas sobre las curvas observadas
#                      'map'    -> Mapas de la predicción y de su varianza, lado
#                                  a lado. Solo si se predijo sobre un mapa.

#' Plot the predictions of a SpatFD object
#'
#' Plots the curves predicted by [predict()][predict.SpatFD] for one variable,
#' either on top of the observed curves or, when the prediction was made over
#' a map, as two maps side by side: the prediction and its variance.
#'
#' @param x An object of class `SpatFDpred`, as returned by
#'   [predict()][predict.SpatFD].
#' @param variable Name of the variable to plot. Defaults to the first
#'   predicted variable.
#' @param mode Either `"curves"`, to draw the predicted curves over the
#'   observed ones, or `"map"`, to draw the maps of the prediction and of its
#'   variance. Defaults to `"map"` when the prediction was made over a map and
#'   to `"curves"` otherwise.
#' @param time Only for `mode = "map"`: time at which the predicted curves are
#'   evaluated to draw the prediction map. Defaults to the first time of the
#'   curves.
#' @param main Title of the plot. By default, for `mode = "map"` it is the
#'   name of the object and the variable, for example `"pred$PM10"`.
#' @param ... Graphical parameters. For `mode = "curves"` they are applied to
#'   the predicted curves (for example `col` or `lwd`); for `mode = "map"`
#'   only `col`, the colour palette of the maps, is used.
#'
#' @details
#' With `mode = "map"` the plot has, from top to bottom: the title; a
#' subtitle in a smaller font with the time, the prediction method and kriging, and whether
#' cross-variograms were used; a colour bar for each map with its scale; the
#' two maps; and under each map what it shows (`Prediction` or `Variance`).
#' The observed sites of the variable are drawn as white points.
#'
#' A predicted curve is a function of time, so a map can only show it at one
#' time: the one given in `time`. The variance map does not depend on `time`,
#' because the prediction variance is that of the whole curve.
#'
#' @return `x`, invisibly. Called for its side effect of drawing a plot.
#'
#' @seealso [predict()][predict.SpatFD]
#'
#' @examples
#' set.seed(1)
#' t <- 1:50
#' dat <- sapply(1:10, function(i) rnorm(1, 5) * sin(2 * pi * t / 50) + rnorm(50, sd = 0.1))
#' colnames(dat) <- paste0("s", 1:10)
#' coords <- cbind(x = runif(10), y = runif(10))
#' rownames(coords) <- colnames(dat)
#' sfd <- spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "sim")
#'
#' sfd <- add_model(sfd, list(
#'   "sim.PC1" = gstat::vgm(psill = 30, model = "Exp", range = 0.5),
#'   "sim.PC2" = gstat::vgm(psill = 0.01, model = "Exp", range = 0.5)
#' ))
#'
#' # Predicted curves over the observed ones
#' pred <- predict(sfd, newcoords = cbind(x = c(0.2, 0.7), y = c(0.5, 0.5)))
#' plot(pred)
#'
#' # Prediction over a map: maps of the prediction and of its variance
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   square <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))))
#'
#'   pred_map <- predict(sfd, newcoords = square, n = 400)
#'   plot(pred_map, time = 25)
#' }
#'
#' @export
plot.SpatFDpred <- function(x, variable = NULL, mode = NULL, time = NULL, main = NULL, ...){

  # Nombre con el que el usuario llamó a la predicción, para el título
  name <- deparse(substitute(x))
  if (length(name) > 1 || nchar(name[1]) > 30) name <- 'prediction'

  if (is.null(variable)) {
    variable <- names(x)[1]
    if (length(x) > 1) message('Plotting the variable ', variable)
  }
  if (!(variable %in% names(x))) {
    stop('Variable ', variable, ' was not predicted. The predicted variables are: ',
         paste(names(x), collapse = ', '), call. = F)
  }

  if (is.null(mode)) mode <- if (is.null(attr(x, 'map'))) 'curves' else 'map'

  # El gráfico no pide <Enter>, aunque el dispositivo esté en modo 'ask'
  old <- par(ask = FALSE)
  on.exit(par(old), add = TRUE)

  switch(mode,
         'curves' = plot_curves(x, variable, main, ...),
         'map' = plot_maps(x, variable, time, main, label = paste0(name, '$', variable), ...),
         stop("'mode' must be the string 'curves' or the string 'map'", call. = F))

  invisible(x)
}

# .............................. FUNCIONES INTERNAS ................................................

# Curvas funcionales ocultas en la predicción
predicted_fd <- function(x, variable) attr(x, 'internal')$pred[[variable]]
observed_fd  <- function(x, variable) attr(x, 'internal')$observed[[variable]]

# GRILLA SOBRE UN MAPA: Cerca de n puntos en una grilla regular dentro del mapa
map_grid <- function(map, n){
  if (!requireNamespace('sf', quietly = TRUE)) {
    stop("The 'sf' package is needed to predict over a map", call. = FALSE)
  }
  stopifnot("'n' must be a single positive number" = is.numeric(n) && length(n) == 1 && n > 0)

  area <- sf::st_union(sf::st_geometry(map))
  grid <- sf::st_coordinates(sf::st_sample(area, size = n, type = 'regular'))[, 1:2, drop = FALSE]
  dimnames(grid) <- NULL
  grid
}

# CURVAS PREDICHAS SOBRE LAS OBSERVADAS: En gris las curvas observadas
# (suavizadas) y en color las predichas. Los parámetros gráficos de ... se
# aplican a las curvas predichas.
plot_curves <- function(x, variable, main = NULL, ...){
  if (is.null(main)) main <- paste(variable, '- Predicted curves')
  style <- utils::modifyList(list(col = 'darkblue', lty = 1), list(...))

  plot(observed_fd(x, variable), col = 'grey80', lty = 1, main = main)
  do.call(lines, c(list(predicted_fd(x, variable)), style))
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
#
# Solo usa las curvas, la varianza y los atributos de la predicción, así que
# no depende del método con que se haya predicho.
plot_maps <- function(x, variable, time, main = NULL, label = variable, ...){
  map <- attr(x, 'map')
  if (is.null(map)) {
    stop("mode = 'map' is only available for predictions made over a map", call. = FALSE)
  }

  curves <- predicted_fd(x, variable)
  times  <- curves$basis$rangeval
  if (is.null(time)) time <- times[1]
  stopifnot(
    "'time' must be a single value within the range of the curves" =
      is.numeric(time) && length(time) == 1 && time >= times[1] && time <= times[2]
  )

  # Título: el que dé el usuario o el nombre de la predicción. Los detalles de
  # cómo se obtuvo van siempre en el subtítulo.
  if (is.null(main)) main <- label
  details <- paste0('Time: ', format(time),
                    ' - Method: ', attr(x, 'method'), ' (', attr(x, 'kriging'), ' kriging)',
                    ' - Crossed: ', if (isTRUE(attr(x, 'crossed'))) 'yes' else 'no')

  grid     <- attr(x, 'newcoords')
  stations <- attr(x, 'internal')$coords[colnames(observed_fd(x, variable)$coefs), , drop = FALSE]
  colours  <- utils::modifyList(list(col = hcl.colors(100, 'viridis', rev = TRUE)), list(...))$col

  # Los puntos forman una grilla regular: se pasan a una matriz para image()
  cell <- function(v) {
    values <- sort(unique(signif(v, 8)))
    if (length(values) > 1) min(diff(values)) else 1
  }
  dx <- cell(grid[, 1])
  dy <- cell(grid[, 2])
  ix <- round((grid[, 1] - min(grid[, 1])) / dx) + 1
  iy <- round((grid[, 2] - min(grid[, 2])) / dy) + 1

  # Cada mapa ocupa dos paneles: su barra de color y, debajo, el mapa
  draw <- function(values, what){
    limits <- range(values)
    if (diff(limits) == 0) limits <- limits + c(-0.5, 0.5)

    # Barra de color en degradado, con la escala de valores
    par(mar = c(2, 6, 0.5, 6))
    image(x = seq(limits[1], limits[2], length.out = length(colours)), y = 1,
          z = matrix(seq_along(colours), ncol = 1), col = colours,
          axes = FALSE, xlab = '', ylab = '')
    axis(1, cex.axis = 1.3, mgp = c(3, 0.6, 0), tcl = -0.3)
    box()

    # Mapa, con el contorno y los sitios observados
    z <- matrix(NA_real_, max(ix), max(iy))
    z[cbind(ix, iy)] <- values

    par(mar = c(3, 1, 0.5, 1))
    image(x = min(grid[, 1]) + dx * (seq_len(max(ix)) - 1),
          y = min(grid[, 2]) + dy * (seq_len(max(iy)) - 1),
          z = z, zlim = limits, col = colours, asp = 1,
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

  draw(as.vector(fda::eval.fd(time, curves)), 'Prediction')
  draw(x[[variable]]$var, 'Variance')
}
