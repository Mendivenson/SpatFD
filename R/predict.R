# ================================== predict() =================================

#' Spatial prediction of functional data
#'
#' Predicts the curves of the variables of a model at unsampled sites, by
#' ordinary kriging (or cokriging, when the model has cross-variograms). The
#' kriging systems were built and checked when the model was created with
#' [spatfd_model()], so predicting only solves them for the new sites.
#'
#' @param object An object of class `SpatFDmodel`, as returned by
#'   [spatfd_model()].
#' @param newcoords Numeric matrix or data frame with the coordinates of the
#'   prediction sites, one row per site. Its columns must be the coordinates
#'   of the variables, matched by name when it has column names and by
#'   position otherwise. Row names, if any, identify the predicted curves.
#'
#'   It can also be a map, an `sf` or `sfc` object of the \pkg{sf} package
#'   with polygons in the same coordinate system as the sites. The curves are
#'   then predicted on a regular grid of points inside the map, see `grid_size`.
#' @param variable Character vector with the names of the variables to
#'   return. Defaults to all the variables of the model. The prediction always
#'   uses every variable of the model; `variable` only selects what is
#'   returned.
#' @param method Prediction method: `"scores"` (the default) or `"lambda"`.
#'   See Details.
#' @param grid_size Only when `newcoords` is a map: approximate number of points of
#'   the grid.
#' @param ... Ignored.
#'
#' @details
#' **Method `"scores"`.** Each score is predicted with its own weights and the
#' curve is rebuilt with the predicted scores:
#' \eqn{\hat\chi_{s_0}(t) = \mu(t) + \sum_k \hat f_k(s_0) \xi_k(t)}.
#'
#' **Method `"lambda"`.** The curve is predicted as a weighted sum of the
#' observed curves, with one weight per observed site and variable, the same
#' for all times. With several variables the curves are combined over the
#' times they share, so `pred` only has those times.
#'
#' Both methods use ordinary kriging: the mean of each field is an unknown
#' constant, so every model can be used, with or without a sill.
#'
#' **Variance.** `var` is the prediction variance of the curve truncated to the
#' components of the model, \eqn{E\|\chi_{s_0} - \hat\chi_{s_0}\|^2} over
#' those components. The variance left out by the truncation is the one not
#' explained by the components, see the element `variance` of the model. A
#' negative variance, which a valid model cannot give, raises a warning.
#'
#' @return An object of class `SpatFDpred`: a list with one element per
#' variable, named after it, each with
#' \describe{
#'   \item{`pred`}{Predicted curves, as a matrix with one row per time of the
#'     observed data and one column per prediction site, the same layout as
#'     the data given to [spatfd()]. The rows are named with the times of the
#'     data (dates, if they were dates).}
#'   \item{`var`}{Prediction variance of the curve at each prediction site.}
#' }
#' The predictions can be plotted with [plot()][plot.SpatFDpred].
#'
#' @section Theory:
#' For a field \eqn{z} observed at \eqn{s_1, \ldots, s_n}, ordinary kriging
#' predicts \eqn{\hat z(s_0) = \sum_i w_i z(s_i)} with weights that add up to
#' one (no bias whatever the mean) and that minimise the prediction error. They
#' solve
#' \deqn{\Gamma w + \mu 1 = \gamma_0, \qquad 1^T w = 1}
#' where \eqn{\Gamma} has the variogram between the observed sites,
#' \eqn{\gamma_0} the variogram between them and \eqn{s_0}, and \eqn{\mu} is a
#' Lagrange multiplier. The prediction variance is
#' \eqn{\gamma_0^T w + \mu}. With several fields (cokriging) the weights of
#' the predicted field add up to one and those of every other field to zero
#' (Bohórquez, Giraldo and Mateu, 2016, 2017).
#'
#' @references
#' Bohórquez, M., Giraldo, R. and Mateu, J. (2016). Optimal sampling for
#' spatial prediction of functional data. *Statistical Methods &
#' Applications*, 25(1), 39-54. \doi{10.1007/s10260-015-0340-9}
#'
#' Bohórquez, M., Giraldo, R. and Mateu, J. (2017). Multivariate functional
#' random fields: prediction and optimal sampling. *Stochastic Environmental
#' Research and Risk Assessment*, 31, 53-70. \doi{10.1007/s00477-016-1266-y}
#'
#' @seealso [spatfd_model()], [plot()][plot.SpatFDpred]
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
#' pred <- predict(model, newcoords = cbind(x = c(0.2, 0.7), y = c(0.5, 0.5)))
#' pred
#' pred$sim$var
#'
#' @export
predict.SpatFDmodel <- function(object,              # Modelo, ver spatfd_model()
                                newcoords,           # Sitios nuevos (matriz, data frame o mapa)
                                variable = NULL,     # Variables que se devuelven
                                method = 'scores',   # Método de predicción
                                grid_size = 5000,    # Puntos de la grilla, solo con un mapa
                                ...){

  stopifnot("'method' must be a single string" = is.character(method) && length(method) == 1)
  method <- match.arg(method, c('scores', 'lambda'))

  data <- get_model_data(object)
  if (is.null(variable)) variable <- names(data)

  unknown <- setdiff(variable, names(data))
  if (length(unknown) > 0) {
    stop('Variables not in the model: ', paste(unknown, collapse = ', '),
         '. The variables of the model are: ', paste(names(data), collapse = ', '), call. = F)
  }

  # MAPAS: En lugar de sitios, una grilla regular de puntos dentro del mapa
  map <- NULL
  if (inherits(newcoords, c('sf', 'sfc'))) {
    map       <- newcoords
    newcoords <- map_grid(map, grid_size)
  }
  newcoords <- check_newcoords(newcoords, colnames(get_coords(data[[1]])))

  # Cada método entrega, por variable, la curva predicha (evaluada y como fd)
  # y su varianza
  curves <- switch(method,
                   'scores' = predict_scores(object, newcoords, variable),
                   'lambda' = predict_lambda(object, newcoords, variable))

  if (any(unlist(lapply(curves, `[[`, 'var')) < 0)) {
    warning('Some prediction variances are negative', call. = F)
  }

  structure(lapply(curves, function(curve) list(pred = curve$pred, var = curve$var)),
            class     = 'SpatFDpred',
            newcoords = newcoords,
            method    = method,
            map       = map,
            internal  = list(pred     = lapply(curves, `[[`, 'fd'),
                             observed = lapply(data[variable], get_fd),
                             coords   = lapply(data[variable], get_coords),
                             time     = lapply(data[variable], get_internal, what = 'time')))
}
