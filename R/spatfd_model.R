# =============================== spatfd_model() ===============================

#' Build a spatial model of functional variables
#'
#' Joins one or several `SpatFD` objects with the variogram models of their
#' scores. The models are compiled once, here, and the result keeps everything
#' that is needed to predict the curves and to validate the model.
#'
#' @param ... Objects of class `SpatFD`, as returned by [spatfd()], with
#'   different names (the attribute `name`). They must be the same objects
#'   whose scores were used to fit `model`.
#' @param model The variogram models of the scores, named with their ids
#'   (`<name>.PC1`, ...; `<id1>.<id2>` for a cross-variogram). One of
#'   \describe{
#'     \item{A fitted `gstat` object}{for example the output of
#'       [gstat::fit.lmc()] on the output of [to_gstat()].}
#'     \item{A named list of `variogramModel` objects}{for example outputs of
#'       [gstat::fit.variogram()], one per score.}
#'     \item{A data frame}{with one row per model structure and the columns
#'       `var` (id), `model` (name of a model in [models]), `sigma` (partial
#'       sill) and, when the model uses them, `rango` (range) and `alpha1`
#'       (shape parameter). The structures of each id are added up.}
#'     \item{A named list of functions}{each one a vectorised variogram with
#'       a single argument, the distance, equal to 0 at distance 0. Its sill
#'       is found from the function itself, see Details.}
#'   }
#'
#' @details
#' **Components.** Each variable is modelled with its first \eqn{K}
#' principal components, where \eqn{K} is given by the direct variograms in
#' `model`: they must be those of `<name>.PC1`, ..., `<name>.PCK`, without gaps.
#' The object keeps a copy of each variable reduced to those \eqn{K}
#' components; `variance` reports the variance they explain.
#'
#' **Cross-variograms** are optional. Two scores without a cross-variogram are
#' treated as independent; with them, they are predicted jointly (cokriging).
#'
#' **Sill.** Models with a sill can be used with simple and ordinary kriging.
#' Models without one (linear with range zero, power and logarithmic) can only
#' be used with ordinary kriging. For a function, the sill is the value of the
#' variogram far from the origin, searched in this order:
#' \enumerate{
#'   \item `f(Inf)`: if finite, it is the sill; if infinite, the model has no
#'     sill.
#'   \item If `f(Inf)` is not defined (for example `sin(Inf)` in a wave
#'     model), its numerical limit at large distances, when the values settle.
#'   \item If the values oscillate without growing (periodic models), their
#'     mean at large distances. If they grow, the model has no sill.
#' }
#' The large distances go from \eqn{10^2} to \eqn{10^{12}} times the size of
#' the area of the sites, so the search does not depend on the units of the
#' coordinates. A message says when the steps 2 or 3 were used.
#'
#' Models of \pkg{gstat} are evaluated with [gstat::variogramLine()], the same
#' code \pkg{gstat} used to fit them. Anisotropic models are not supported yet.
#' The function does not check that `model` was fitted on the scores of the
#' objects given.
#'
#' @return An object of class `SpatFDmodel`: a list with
#' \describe{
#'   \item{`variance`}{Data frame with one row per variable: basis type,
#'     number of components in the model and the proportion of variance they
#'     explain.}
#'   \item{`variograms`}{Data frame with one row per variogram: id, type
#'     (direct or cross), where it came from, its structures and its sill
#'     (`NA` without a sill).}
#' }
#'
#' @section Theory:
#' Prediction needs the covariances of the scores, which come from their
#' variograms:
#' \deqn{C(h) = C(0) - \gamma(h)}
#' where \eqn{C(0)} is the sill. Without a sill, ordinary kriging can still be
#' used, because its constraints cancel any constant added to the covariance
#' (Chilès and Delfiner, 2012). The model describes the process truncated to
#' its first \eqn{K} components per variable: the variance left out is the one
#' not explained by them (Bohórquez, Giraldo and Mateu, 2016, 2017).
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
#' Chilès, J.-P. and Delfiner, P. (2012). *Geostatistics: Modeling Spatial
#' Uncertainty* (2nd ed.). John Wiley & Sons. \doi{10.1002/9781118136188}
#'
#' @seealso [to_gstat()], [models], [plot()][plot.SpatFDmodel]
#'
#' @examples
#' set.seed(1)
#' dat <- sapply(1:10, function(i) rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1))
#' colnames(dat) <- paste0("s", 1:10)
#' coords <- cbind(x = runif(10), y = runif(10))
#' rownames(coords) <- colnames(dat)
#' sim <- spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "sim")
#'
#' # A data frame of models from the catalog
#' table <- data.frame(var = c("sim.PC1", "sim.PC1"), model = c("nugget", "exponential"),
#'                     rango = c(NA, 0.5), sigma = c(0.1, 2))
#' m <- spatfd_model(sim, model = table)
#' m
#'
#' @examplesIf requireNamespace("gstat", quietly = TRUE)
#' # The same model as a gstat variogram model
#' summary(spatfd_model(sim, model = list(sim.PC1 = gstat::vgm(2, "Exp", 0.5, nugget = 0.1))))
#'
#' @export
spatfd_model <- function(...,                # Objetos SpatFD
                         model){             # Modelos de los variogramas de los scores

  if (missing(model)) stop("The variogram models must be given in 'model'", call. = F)

  objects  <- spatfd_list(...)

  # Tamaño del área de los sitios: la diagonal del rectángulo que los contiene
  coords   <- do.call(rbind, lapply(objects, get_coords))
  scale    <- sqrt(sum(apply(coords, 2, function(v) diff(range(v)))^2))

  compiled <- normalize_model(model, scale)           # Una entrada por id, ver spatfd_model_UTILS.R
  ids      <- match_ids(names(compiled), objects)
  K        <- deduce_k(ids, objects)

  # COPIAS: Cada variable con las componentes que usa el modelo
  data     <- mapply(truncate_spatfd, objects, K, SIMPLIFY = F)
  variance <- do.call(rbind, lapply(data, `[[`, 'variance'))

  variograms <- data.frame('variogram'  = ids$id,
                           'type'       = ids$type,
                           'source'     = vapply(compiled, `[[`, character(1), 'source'),
                           'structures' = vapply(compiled, `[[`, character(1), 'structures'),
                           'sill'       = vapply(compiled, `[[`, numeric(1), 'sill'),
                           row.names = NULL)

  # SISTEMAS DE KRIGING: Se arman, validan y factorizan una sola vez, para los
  # dos métodos de predicción. Un modelo no válido se detiene aquí.
  check_coregionalization(compiled, ids, data)
  functions <- lapply(compiled, `[[`, 'variogram')

  structure(list(variance = variance, variograms = variograms),
            class    = 'SpatFDmodel',
            data     = data,
            internal = list(input      = model,         # Para rearmar el modelo sin un sitio, ver leave_one_out()
                            ids        = ids,
                            variograms = functions,
                            sills      = vapply(compiled, `[[`, numeric(1), 'sill'),
                            scores     = scores_systems(data, ids, functions),
                            lambda     = lambda_systems(data, ids, functions)))
}
