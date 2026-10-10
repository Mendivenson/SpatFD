# ================================ predict() ===================================

#' Spatial prediction of functional data
#'
#' Predicts the curves of one or several functional variables at unsampled
#' sites. With the default method, the scores of the principal components are
#' predicted by simple kriging (or simple cokriging) and the curves are rebuilt
#' from them; with the method `"lambda"`, the curves are predicted as a linear
#' combination of the observed curves. The systems are solved inside the
#' package, in matrix form.
#'
#' @param object An object of class `SpatFD`, as returned by [spatfd()].
#' @param newcoords Numeric matrix or data frame with the coordinates of the
#'   prediction sites, one row per site. Its columns must be the coordinates
#'   of `object`, matched by name when it has column names and by position
#'   otherwise. Row names, if any, identify the predicted curves.
#'
#'   It can also be a map, an `sf` or `sfc` object of the \pkg{sf} package
#'   with polygons in the same coordinate system as the sites of `object`. The
#'   curves are then predicted on a regular grid of points inside the map, see
#'   `n`.
#' @param model Name of one of the models registered in `object` with
#'   [add_model()]; `object$models` lists them. It can be omitted when
#'   the object has a single model. The model must have a variogram for every
#'   score of the variables in `variable`; cross-variograms between two scores
#'   are optional, see Details.
#' @param variable Character vector with the names of the variables to
#'   predict. Defaults to all the variables covered by the model. With
#'   several variables their scores are predicted jointly.
#' @param method Prediction method, a single one per call: `"scores"` (the
#'   default) or `"lambda"`. See Details and Theory.
#' @param n Only when `newcoords` is a map: approximate number of points of
#'   the grid on which the curves are predicted.
#' @param plot Only when `newcoords` is a map: if `TRUE` (the default), the
#'   maps of the prediction and of its variance are drawn side by side, see
#'   [plot.SpatFDpred()].
#' @param ... Ignored.
#'
#' @details
#' The model defines the system that is solved. Two scores are predicted
#' jointly only when the model has a cross-variogram for them; otherwise they
#' are treated as independent. So a model with only direct variograms gives
#' one kriging per score, and a model with cross-variograms, such as a linear
#' model of coregionalisation fitted with [gstat::fit.lmc()], gives cokriging.
#'
#' The models are compiled once, when they are registered, and reused in
#' every prediction. If the model cannot be used, the error shows the summary
#' table of the models of the object.
#'
#' Simple kriging works with covariances, so every model must have a sill:
#' `linear`, `power`, `logarithmic`, `spline` and `legendre` cannot be used.
#'
#' The method `"lambda"` reproduces the predictor of the CRAN version of the
#' package and has, for now, these limitations:
#' * It predicts a single variable with a model without cross-variograms; any
#'   other case is an error.
#' * The weights are applied to the observed curves as they are, without
#'   removing their mean first.
#' * The variance of each score is taken from the eigenvalue of its principal
#'   component, not from the sill of the fitted model, so the prediction
#'   variance can be negative when the two differ much.
#'
#' @return An object of class `SpatFDpred`: a list with one element per
#' predicted variable, named after it, each with
#' \describe{
#'   \item{`pred`}{Predicted curves, as a functional data object of \pkg{fda}
#'     with one curve per prediction site.}
#'   \item{`var`}{Prediction variance of the curve at each prediction site:
#'     the sum of the prediction variances of its scores.}
#' }
#' The list always has this form, also when a single variable is predicted.
#' The coordinates of the prediction sites are kept in the attribute
#' `newcoords`. When the prediction is made over a map and the maps are
#' drawn, the object is returned invisibly.
#'
#' The predictions can be plotted with [plot()][plot.SpatFDpred].
#'
#' @section Theory:
#' The scores are centred, so they are predicted by simple kriging with a
#' known mean of zero. Let \eqn{f} be the vector with the observed scores, of
#' all the components and variables involved, \eqn{\Sigma} its covariance
#' matrix and \eqn{c_0} its covariance with the score \eqn{f_k(s_0)} at the
#' unsampled site \eqn{s_0}. The predictor and its variance are
#' \deqn{\hat{f}_k(s_0) = c_0^T \Sigma^{-1} f, \qquad
#'   \sigma_k^2(s_0) = C_k(0) - c_0^T \Sigma^{-1} c_0}
#' The covariances come from the fitted variograms, \eqn{C(h) = C(0) -
#' \gamma(h)}, where \eqn{C(0)} is the sum of the partial sills.
#'
#' \eqn{\Sigma} is made of one block per pair of scores. Without
#' cross-variograms it is block diagonal and the system splits into one
#' kriging per score; with them, it is the simple cokriging system. The curve
#' is then rebuilt as (Bohorquez, Giraldo and Mateu, 2016, 2017)
#' \deqn{\hat{\chi}_{s_0}(t) = \mu(t) + \sum_{k=1}^{K} \hat{f}_k(s_0) \xi_k(t)}
#' Because the eigenfunctions are orthonormal, the sum of the prediction
#' variances of the scores measures the uncertainty of the whole curve.
#'
#' The method `"lambda"` predicts the curve directly as a linear combination
#' of the observed curves, with one scalar weight per observed site (Giraldo,
#' Delicado and Mateu, 2011):
#' \deqn{\hat{\chi}_{s_0}(t) = \sum_{i=1}^{n} \lambda_i \chi_{s_i}(t)}
#' The weights minimise \eqn{E \| \chi_{s_0} - \hat{\chi}_{s_0} \|^2} and
#' solve \eqn{\Omega \lambda = c}, where \eqn{\Omega = \sum_k \Omega_k} and
#' \eqn{c = \sum_k c_k} add up, over the scores, the covariance matrices
#' between observed sites and the covariances with \eqn{s_0}. The prediction
#' variance is
#' \deqn{\sum_{k=1}^{K} \left( \eta_k - 2 c_k^T \lambda +
#'   \lambda^T \Omega_k \lambda \right)}
#' with \eqn{\eta_k} the eigenvalue of the \eqn{k}-th principal component.
#' Unlike the method `"scores"`, the same weights are used for all the scores
#' and the curves are not truncated to the first \eqn{K} components.
#'
#' @references
#' Bohorquez, M., Giraldo, R. and Mateu, J. (2016). Optimal sampling for
#' spatial prediction of functional data. *Statistical Methods &
#' Applications*, 25(1), 39-54. \doi{10.1007/s10260-015-0340-9}
#'
#' Bohorquez, M., Giraldo, R. and Mateu, J. (2017). Multivariate functional
#' random fields: prediction and optimal sampling. *Stochastic Environmental
#' Research and Risk Assessment*, 31, 53-70. \doi{10.1007/s00477-016-1266-y}
#'
#' Chilès, J.-P. and Delfiner, P. (2012). *Geostatistics: Modeling Spatial
#' Uncertainty* (2nd ed.). John Wiley & Sons. \doi{10.1002/9781118136188}
#'
#' Giraldo, R., Delicado, P. and Mateu, J. (2011). Ordinary kriging for
#' function-valued spatial data. *Environmental and Ecological Statistics*,
#' 18, 411-426. \doi{10.1007/s10651-010-0143-y}
#'
#' @seealso [add_model()] to register models in the object, and
#'   [plot.SpatFDpred()] to plot the predictions.
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
#' model <- toSpatFD(list(
#'   "sim.PC1" = gstat::vgm(psill = 30, model = "Exp", range = 0.5),
#'   "sim.PC2" = gstat::vgm(psill = 0.01, model = "Exp", range = 0.5)
#' ))
#'
#' add_model(sfd, model, name = "exponential")
#' sfd$models
#'
#' pred <- predict(sfd, newcoords = cbind(x = c(0.2, 0.7), y = c(0.5, 0.5)), model = "exponential")
#' pred
#' pred$sim$var
#' plot(pred)
#'
#' # The same prediction with the method lambda
#' pred_lambda <- predict(sfd, newcoords = cbind(x = c(0.2, 0.7), y = c(0.5, 0.5)),
#'                        model = "exponential", method = "lambda")
#' plot(pred_lambda)
#'
#' @importFrom utils capture.output
#' @export
predict.SpatFD <- function(object, newcoords, model = NULL, variable = NULL, method = 'scores',
                           n = 5000, plot = TRUE, ...){

  stopifnot("'method' must be a single string" = is.character(method) && length(method) == 1)
  method <- match.arg(method, c('scores', 'lambda'))

  # MAPAS: Si en lugar de sitios se recibe un mapa, se predice sobre una grilla
  # regular de puntos dentro de él
  map <- NULL
  if (inherits(newcoords, c('sf', 'sfc'))) {
    map       <- newcoords
    newcoords <- map_grid(map, n)
  }

  # Cada método entrega las curvas predichas y su varianza; los pasos
  # intermedios (scores predichos, pesos) no se le entregan al usuario
  result <- switch(method,
                   'scores' = predict_scores(object, newcoords, model, variable),
                   'lambda' = predict_lambda(object, newcoords, model, variable))

  # Además de las predicciones se guarda lo necesario para graficarlas: el mapa
  # (si lo hay) y, por referencia, las curvas y los sitios observados
  # También cómo se obtuvieron: el método y si se usaron variogramas cruzados
  result <- structure(result$curves,
                      class = 'SpatFDpred',
                      newcoords = result$newcoords,
                      method = method,
                      crossed = result$crossed,
                      map = map,
                      fd = attr(object, 'fd'),
                      coords = attr(object, 'data')$coords)

  if (!is.null(map) && plot) {
    # Aquí la predicción todavía no tiene nombre: el título lleva el de la variable
    if (length(result) > 1) message('Plotting the variable ', names(result)[1])
    plot_maps(result, names(result)[1], time = NULL, main = NULL, label = names(result)[1])
    return(invisible(result))
  }
  result
}

# FUNCIÓN DE IMPRESIÓN: Las curvas predichas (fd) no se imprimen, solo la
#                       varianza de predicción en los primeros sitios

#' @rdname predict.SpatFD
#' @param x An object of class `SpatFDpred`, as returned by `predict()`.
#' @export
print.SpatFDpred <- function(x, ...){
  sites <- nrow(attr(x, 'newcoords'))
  shown <- min(sites, 10)

  cat('\n\n-------------------------------------------------\n')
  cat('Predicted curves at ', sites, ' sites', if (!is.null(attr(x, 'map'))) ' over a map',
      ' for: ', paste(names(x), collapse = ', '), '\n', sep = '')
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
#' subtitle in a smaller font with the time, the prediction method and whether
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
#' add_model(sfd, list(
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
#' @importFrom graphics axis box image layout lines mtext plot.new points text
#' @importFrom grDevices hcl.colors
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

  switch(mode,
         'curves' = plot_curves(x, variable, main, ...),
         'map' = plot_maps(x, variable, time, main, label = paste0(name, '$', variable), ...),
         stop("'mode' must be the string 'curves' or the string 'map'", call. = F))

  invisible(x)
}

# .............................. FUNCIONES INTERNAS ................................................

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

  plot(attr(x, 'fd')[[variable]], col = 'grey80', lty = 1, main = main)
  do.call(lines, c(list(x[[variable]]$pred), style))
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

  curves <- x[[variable]]$pred
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
                    ' - Method: ', attr(x, 'method'),
                    ' - Crossed: ', if (isTRUE(attr(x, 'crossed'))) 'yes' else 'no')

  grid     <- attr(x, 'newcoords')
  stations <- attr(x, 'coords')[colnames(attr(x, 'fd')[[variable]]$coefs), , drop = FALSE]
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
    graphics::par(mar = c(2, 6, 0.5, 6))
    image(x = seq(limits[1], limits[2], length.out = length(colours)), y = 1,
          z = matrix(seq_along(colours), ncol = 1), col = colours,
          axes = FALSE, xlab = '', ylab = '')
    axis(1, cex.axis = 1.3, mgp = c(3, 0.6, 0), tcl = -0.3)
    box()

    # Mapa, con el contorno y los sitios observados
    z <- matrix(NA_real_, max(ix), max(iy))
    z[cbind(ix, iy)] <- values

    graphics::par(mar = c(3, 1, 0.5, 1))
    image(x = min(grid[, 1]) + dx * (seq_len(max(ix)) - 1),
          y = min(grid[, 2]) + dy * (seq_len(max(iy)) - 1),
          z = z, zlim = limits, col = colours, asp = 1,
          axes = FALSE, xlab = '', ylab = '')
    plot(sf::st_geometry(map), add = TRUE)
    points(stations, pch = 21, bg = 'white', cex = 1.2)

    mtext(what, side = 1, line = 1.2, cex = 1.3, font = 2)
  }

  old <- graphics::par(mar = c(0, 0, 0, 0))
  on.exit({ layout(1); graphics::par(old) })

  layout(matrix(c(1, 1,
                  2, 4,
                  3, 5), nrow = 3, byrow = TRUE), heights = c(2, 1.3, 10))

  plot.new()
  text(0.5, 0.68, main, cex = 2.4, font = 2)
  text(0.5, 0.18, details, cex = 1.5)

  draw(as.vector(fda::eval.fd(time, curves)), 'Prediction')
  draw(x[[variable]]$var, 'Variance')
}

# MÉTODO SCORES: Kriging (o cokriging) de los scores y reconstrucción de las
# curvas como media + combinación lineal de las funciones propias, todo sobre
# los coeficientes de la base.
#
# Los dos métodos devuelven lo mismo, que es lo único que usa predict():
#   curves    -> Lista por variable con las curvas (pred) y su varianza (var)
#   newcoords -> Coordenadas de los sitios nuevos, ya verificadas
#   crossed   -> Si se usó algún variograma cruzado
predict_scores <- function(object, newcoords, model, variable){
  kriging <- krige_scores(object, newcoords, model, variable)
  sites   <- rownames(kriging$newcoords)

  curves <- lapply(names(kriging$pcs), function(v) {
    pca <- attr(object, 'pca')[[v]]
    pcs <- kriging$pcs[[v]]

    coefs <- pca$harmonics$coefs %*% t(kriging$pred[, pcs, drop = FALSE])
    coefs <- coefs + as.vector(pca$meanfd$coefs)
    colnames(coefs) <- sites

    fdnames      <- attr(object, 'fd')[[v]]$fdnames
    fdnames[[2]] <- sites

    list(pred = fda::fd(coefs, pca$harmonics$basis, fdnames),
         var  = rowSums(kriging$variance[, pcs, drop = FALSE]))
  })
  names(curves) <- names(kriging$pcs)

  list(curves = curves, newcoords = kriging$newcoords, crossed = kriging$crossed)
}

# MÉTODO LAMBDA: La curva predicha es una combinación lineal de las curvas
# observadas, con un único juego de pesos por sitio nuevo. Los pesos salen de
# la covarianza de las curvas completas, que es la suma de las covarianzas de
# los scores. Por ahora solo para una variable y sin variogramas cruzados.
predict_lambda <- function(object, newcoords, model, variable){
  setup <- prediction_setup(object, newcoords, model, variable)

  crossed <- !all(setup$compiled$table$var %in% names(score_ids(object)))
  if (length(setup$pcs) > 1 || crossed) {
    stop("The method 'lambda' only predicts one variable with a model without cross-variograms. ",
         if (crossed) paste0("The model '", setup$model, "' has cross-variograms")
         else paste0("Choose one variable with 'variable': ", paste(names(setup$pcs), collapse = ', ')),
         '\nThe models of the object are:\n', setup$summary_table(), call. = FALSE)
  }

  v      <- names(setup$pcs)
  ids    <- setup$pcs[[v]]
  places <- setup$sites[[1]]              # Los scores de una variable comparten sitios
  covs   <- score_covariances(setup$compiled, ids)

  # Covarianza de cada score entre sitios observados (Omega) y con los sitios
  # nuevos (C0). Las de las curvas completas son la suma sobre los scores.
  d  <- cross_dist(places, places)
  d0 <- cross_dist(places, setup$newcoords)

  evaluate <- function(covariance, d) matrix(covariance(d), nrow(d), ncol(d))
  Omega <- lapply(ids, function(id) evaluate(covs[[id, id]], d))
  C0    <- lapply(ids, function(id) evaluate(covs[[id, id]], d0))

  R <- tryCatch(chol(Reduce(`+`, Omega)), error = function(e) {
    stop('The covariance matrix of the curves is not positive definite. ',
         'Check for duplicated sites or for an invalid model', call. = FALSE)
  })
  weights <- backsolve(R, backsolve(R, Reduce(`+`, C0), transpose = TRUE))   # Sitios x sitios nuevos

  # VARIANZA: Suma sobre los scores de  Var(f_k) - 2 c_k' lambda + lambda' Omega_k lambda.
  # Como en la versión CRAN, Var(f_k) es el valor propio del componente y no
  # la meseta del modelo ajustado, por lo que puede salir negativa.
  eigenvalues <- attr(object, 'pca')[[v]]$values[seq_along(ids)]
  variance <- Reduce(`+`, lapply(seq_along(ids), function(k) {
    eigenvalues[k] - 2 * colSums(C0[[k]] * weights) + colSums(weights * (Omega[[k]] %*% weights))
  }))
  names(variance) <- rownames(setup$newcoords)

  # CURVAS: Como en la versión CRAN, los pesos se aplican a las curvas
  # observadas sin centrar. Los pesos son de kriging simple y no suman uno, así
  # que una implementación más cuidadosa los aplicaría a las curvas centradas y
  # sumaría la media al final; basta con cambiar las dos líneas marcadas.
  observed <- attr(object, 'fd')[[v]]
  coefs    <- observed$coefs[, rownames(places), drop = FALSE]    # <- restar aquí la media
  coefs    <- coefs %*% weights                                   # <- sumar aquí la media
  colnames(coefs) <- rownames(setup$newcoords)

  fdnames      <- observed$fdnames
  fdnames[[2]] <- rownames(setup$newcoords)

  curves <- list(list(pred = fda::fd(coefs, observed$basis, fdnames), var = variance))
  names(curves) <- v

  list(curves = curves, newcoords = setup$newcoords, crossed = FALSE)
}

# PREPARACIÓN COMÚN A LOS DOS MÉTODOS: Verifica los argumentos y deja listo lo
# que se necesita para predecir:
#   model, compiled -> Nombre del modelo escogido y su versión compilada
#   pcs             -> Lista con los scores que pertenecen a cada variable
#   sites, values   -> Sitios y valores observados de cada score
#   newcoords       -> Coordenadas de los sitios nuevos, ya verificadas
#   summary_table   -> Función que devuelve la tabla de modelos, para los errores
prediction_setup <- function(object, newcoords, model = NULL, variable = NULL){

  # ============================ MODELO A USAR =================================
  # Los modelos ya están compilados dentro del objeto, ver add_model(). En
  # caso de error se muestra la tabla resumen de los modelos registrados.
  available     <- ls(attr(object, 'models'))
  summary_table <- function() paste(capture.output(print(object$models)), collapse = '\n')

  if (length(available) == 0) {
    stop('The object has no models. Register one with add_model()', call. = FALSE)
  }
  if (is.null(model)) {
    if (length(available) > 1) {
      stop("The object has several models, choose one with 'model':\n", summary_table(), call. = FALSE)
    }
    model <- available
  }
  if (!is.character(model) || length(model) != 1 || !(model %in% available)) {
    stop("'model' must be the name of one of the models of the object:\n", summary_table(),
         call. = FALSE)
  }
  compiled <- attr(object, 'models')[[model]]

  # Scores del objeto (un campo = un score) y variable a la que pertenecen
  all_ids <- score_ids(object)

  if (is.null(variable)) {
    variable <- unique(all_ids[intersect(names(all_ids), names(compiled$variograms))])
  }

  coords <- attr(object, 'data')$coords

  stopifnot(
    'Variables not found in object' = all(variable %in% var_names(object)),
    "'newcoords' must be a numeric matrix or data frame without missing values" =
      (is.matrix(newcoords) || is.data.frame(newcoords)) &&
      is.numeric(as.matrix(newcoords)) && !anyNA(newcoords) && nrow(newcoords) > 0
  )

  # Coordenadas de predicción: por nombre si los tiene, por posición si no
  newcoords <- as.matrix(newcoords)
  if (all(colnames(coords) %in% colnames(newcoords))) {
    newcoords <- newcoords[, colnames(coords), drop = FALSE]
  } else if (is.null(colnames(newcoords)) && ncol(newcoords) == ncol(coords)) {
    colnames(newcoords) <- colnames(coords)
  } else {
    stop("'newcoords' must have the coordinates of the object: ",
         paste(colnames(coords), collapse = ', '), call. = FALSE)
  }
  if (is.null(rownames(newcoords))) {
    rownames(newcoords) <- paste0('new', seq_len(nrow(newcoords)))
  }

  # ============================ SCORES A PREDECIR =============================
  pcs <- lapply(variable, function(v) names(all_ids)[all_ids == v])
  names(pcs) <- variable
  ids <- unlist(pcs, use.names = FALSE)

  missing_models <- setdiff(ids, names(compiled$variograms))
  if (length(missing_models) > 0) {
    stop("The model '", model, "' has no variogram for the scores: ",
         paste(missing_models, collapse = ', '),
         '\nThe models of the object are:\n', summary_table(), call. = FALSE)
  }

  # Sitios con score faltante se descartan campo a campo: las variables pueden
  # estar observadas en conjuntos de sitios distintos
  scores   <- object$scores
  observed <- lapply(ids, function(id) !is.na(scores[, id]))
  sites    <- lapply(observed, function(obs) scores[obs, colnames(coords), drop = FALSE])
  values   <- Map(function(id, obs) scores[obs, id], ids, observed)

  list(model = model, compiled = compiled, pcs = pcs, sites = sites, values = values,
       newcoords = newcoords, summary_table = summary_table)
}

# KRIGING DE LOS SCORES: Predice por kriging simple los scores de las variables
# pedidas en los sitios nuevos. Es el paso común a todo lo que necesite los
# scores predichos (método scores, validación cruzada, ...) y devuelve:
#   pred      -> Matriz de scores predichos (sitios nuevos x scores)
#   variance  -> Matriz con la varianza de predicción de cada score
#   pcs       -> Lista con los scores que pertenecen a cada variable
#   newcoords -> Coordenadas de los sitios nuevos, ya verificadas
#   crossed   -> Si se usó algún variograma cruzado entre los scores predichos
krige_scores <- function(object, newcoords, model = NULL, variable = NULL){
  setup <- prediction_setup(object, newcoords, model, variable)
  ids   <- unlist(setup$pcs, use.names = FALSE)

  covs    <- score_covariances(setup$compiled, ids)
  kriging <- simple_kriging(setup$sites, setup$values, covs, setup$newcoords)
  dimnames(kriging$pred) <- dimnames(kriging$variance) <- list(rownames(setup$newcoords), ids)

  # ¿Se usó algún variograma cruzado entre los scores predichos?
  crossed <- !all(vapply(covs[upper.tri(covs)], is.null, logical(1)))

  list(pred = kriging$pred, variance = kriging$variance, pcs = setup$pcs,
       newcoords = setup$newcoords, crossed = crossed)
}

# COVARIANZAS DE LOS SCORES: A partir de un modelo ya compilado, ver
# add_model(), se construye una matriz-lista con la función de covarianza
# de cada pareja de campos, C(h) = meseta - gamma(h). Las parejas sin
# variograma cruzado quedan en NULL, es decir, se toman como independientes.
score_covariances <- function(compiled, ids){
  unbounded <- c('linear', 'power', 'logarithmic', 'spline', 'legendre')

  q    <- length(ids)
  covs <- matrix(vector('list', q * q), q, q, dimnames = list(ids, ids))

  for (i in seq_len(q)) {
    for (j in i:q) {
      # El id de un variograma cruzado puede venir en cualquiera de los dos órdenes
      key <- if (i == j) ids[i] else paste(ids[c(i, j)], ids[c(j, i)], sep = '.')
      key <- intersect(key, names(compiled$variograms))
      if (length(key) == 0) next

      structures <- compiled$table$model[compiled$table$var == key[1]]
      if (any(structures %in% unbounded)) {
        stop("The model of '", key[1], "' has no sill (",
             paste(intersect(structures, unbounded), collapse = ', '),
             '), so it cannot be used for simple kriging', call. = FALSE)
      }
      if (is.na(compiled$sills[[key[1]]])) {
        stop("The model of '", key[1], "' has missing values in 'sigma'", call. = FALSE)
      }

      covs[[i, j]] <- covs[[j, i]] <- local({
        sill  <- compiled$sills[[key[1]]]
        gamma <- compiled$variograms[[key[1]]]
        function(h) sill - gamma(h)
      })
    }
  }
  covs
}

# DISTANCIAS ENTRE DOS CONJUNTOS DE SITIOS: Se calcula coordenada a coordenada
# para que la distancia entre dos sitios iguales sea exactamente cero (el
# efecto pepita depende de ello).
cross_dist <- function(a, b){
  d2 <- matrix(0, nrow(a), nrow(b))
  for (k in seq_len(ncol(a))) {
    d2 <- d2 + outer(a[, k], b[, k], '-')^2
  }
  sqrt(d2)
}

# SOLVER DE KRIGING SIMPLE (media conocida e igual a cero):
#   sites     -> Lista con las coordenadas de los sitios de cada campo
#   values    -> Lista con los valores observados de cada campo
#   covs      -> Matriz-lista de funciones de covarianza, ver score_covariances()
#   newcoords -> Coordenadas de los sitios de predicción
#
# Sigma se arma por bloques, uno por pareja de campos. Sin covarianzas cruzadas
# es diagonal por bloques y el sistema equivale a un kriging por campo; con
# ellas es el sistema de cokriging. Los casos univariado y multivariado son
# entonces el mismo sistema.
#
# Los campos que no están conectados por ninguna covarianza cruzada no
# comparten información, así que cada grupo de campos conectados se resuelve
# por aparte: da lo mismo que resolver el sistema completo y evita factorizar
# una matriz llena de ceros.
simple_kriging <- function(sites, values, covs, newcoords){
  q <- length(sites)
  m <- nrow(newcoords)

  # Los modelos no siempre conservan las dimensiones de la matriz de distancias.
  # Si las distancias son entre un conjunto de sitios y él mismo la matriz es
  # simétrica y basta evaluar la mitad.
  evaluate <- function(covariance, d, symmetric = FALSE){
    if (!symmetric) return(matrix(covariance(d), nrow(d), ncol(d)))

    half <- lower.tri(d, diag = TRUE)
    out  <- matrix(0, nrow(d), ncol(d))
    out[half] <- covariance(d[half])
    out[upper.tri(out)] <- t(out)[upper.tri(out)]
    out
  }

  # DISTANCIAS: Los scores de una misma variable comparten sitios, así que cada
  # matriz de distancias se calcula una sola vez por conjunto distinto de sitios
  places <- list()
  place  <- integer(q)                   # Conjunto de sitios de cada campo
  for (i in seq_len(q)) {
    known <- Position(function(p) identical(p, sites[[i]]), places)
    if (is.na(known)) {
      places <- c(places, list(sites[[i]]))
      known  <- length(places)
    }
    place[i] <- known
  }

  computed <- list()
  distance <- function(a, b){            # b = 0 son los sitios de predicción
    key <- paste(a, b)
    if (is.null(computed[[key]])) {
      computed[[key]] <<- cross_dist(places[[a]], if (b == 0) newcoords else places[[b]])
    }
    computed[[key]]
  }

  # GRUPOS: Campos conectados, directa o indirectamente, por covarianzas cruzadas
  group <- seq_len(q)
  for (i in seq_len(q)) {
    for (j in seq_len(q)) {
      if (!is.null(covs[[i, j]])) group[group == group[j]] <- group[i]
    }
  }

  pred <- variance <- matrix(0, m, q)

  for (g in unique(group)) {
    fields <- which(group == g)
    k <- length(fields)
    n <- vapply(sites[fields], nrow, integer(1))

    rows <- split(seq_len(sum(n)), rep(seq_len(k), n))         # Filas de cada campo
    cols <- split(seq_len(m * k), rep(seq_len(k), each = m))   # Columnas de cada campo objetivo

    Sigma <- matrix(0, sum(n), sum(n))    # Covarianza entre observaciones
    C0    <- matrix(0, sum(n), m * k)     # Covarianza entre observaciones y sitios nuevos

    for (a in seq_len(k)) {
      for (b in a:k) {
        i <- fields[a]
        j <- fields[b]
        if (is.null(covs[[i, j]])) next
        same_sites <- place[i] == place[j]

        # Sigma es simétrica: solo se evalúa un bloque de cada pareja de campos
        block <- evaluate(covs[[i, j]], distance(place[i], place[j]), symmetric = same_sites)
        Sigma[rows[[a]], rows[[b]]] <- block
        if (a < b) Sigma[rows[[b]], rows[[a]]] <- t(block)

        # La covarianza con los sitios nuevos solo depende de los sitios del
        # campo observado: si ambos campos comparten sitios, es la misma
        block <- evaluate(covs[[i, j]], distance(place[i], 0))
        C0[rows[[a]], cols[[b]]] <- block
        if (a < b) {
          C0[rows[[b]], cols[[a]]] <- if (same_sites) block else evaluate(covs[[i, j]], distance(place[j], 0))
        }
      }
    }

    # Una sola factorización sirve para todos los sitios y campos del grupo.
    # Cholesky falla si Sigma no es definida positiva, lo que delata un modelo
    # inválido.
    R <- tryCatch(chol(Sigma), error = function(e) {
      stop('The covariance matrix of the scores is not positive definite. ',
           'Check for duplicated sites or for an invalid model', call. = FALSE)
    })
    #
    # Con Sigma = R'R y V = R'^-1 C0 no hace falta calcular los pesos Sigma^-1 C0:
    #   predicción = C0' Sigma^-1 z  = V' (R'^-1 z)
    #   varianza   = C(0) - C0' Sigma^-1 C0 = C(0) - suma por columnas de V^2
    V <- backsolve(R, C0, transpose = TRUE)
    z <- backsolve(R, unlist(values[fields]), transpose = TRUE)

    sills <- vapply(fields, function(j) covs[[j, j]](0), numeric(1))

    pred[, fields]     <- matrix(crossprod(V, z), m, k)
    variance[, fields] <- matrix(rep(sills, each = m) - colSums(V^2), m, k)
  }

  list(pred = pred, variance = variance)
}
