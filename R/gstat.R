# ======================= INTERCAMBIO CON gstat ================================
# El ajuste de los variogramas se hace por fuera del paquete. togstat() exporta
# los scores y toSpatFD() trae de vuelta los modelos ajustados.

#' Export the scores of a SpatFD object to gstat
#'
#' Builds a `gstat` object with one entry per principal component score of the
#' selected variables, each with a constant mean (`score ~ 1`). The result can
#' be used with [gstat::variogram()], [gstat::fit.variogram()] or
#' [gstat::fit.lmc()].
#'
#' @param sfd An object of class `SpatFD`, as returned by [spatfd()].
#' @param vars Character vector with the names of the variables to export.
#'   Defaults to all the variables in `sfd`.
#'
#' @details
#' Sites with a missing score are dropped from the corresponding entry, so
#' variables observed at different sets of sites can be exported together.
#'
#' @return An object of class `gstat` whose ids are the names of the score
#'   columns (`<variable>.PC1`, `<variable>.PC2`, ...).
#'
#' @section Theory:
#' After the functional principal component analysis, the spatial dependence
#' between curves is determined by the spatial dependence of the scores
#' \eqn{f_k(s)}, which are scalar random fields. Each score is therefore
#' treated as a classical geostatistical variable: one variogram per score in
#' the univariate case, plus the cross-variograms between scores when several
#' functional variables are modelled together (Bohorquez, Giraldo and Mateu,
#' 2016, 2017). \pkg{gstat} (Pebesma, 2004) is one tool to estimate and fit
#' them.
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
#' Pebesma, E. J. (2004). Multivariable geostatistics in S: the gstat package.
#' *Computers & Geosciences*, 30, 683-691. \doi{10.1016/j.cageo.2004.03.012}
#'
#' @seealso [toSpatFD()] to translate the fitted models back.
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
#' g <- togstat(sfd)
#' names(g$data)
#'
#' @export
togstat <- function(sfd, vars = var_names(sfd)) {

  stopifnot(
    'function only available for SpatFD objects' = inherits(sfd, 'SpatFD'),
    'Variables not found in object' = all(vars %in% var_names(sfd))
  )

  # Scores de las variables escogidas, con las coordenadas de sus sitios. Se
  # conserva el orden de sfd$scores: el ajuste de gstat::fit.lmc() depende de él
  ids    <- score_ids(sfd)
  keep   <- colnames(sfd$scores) %in% names(ids)[ids %in% vars]
  scores <- as.data.frame(sfd$scores[, keep, drop = FALSE])
  sp::coordinates(scores) <- get_coords(sfd)[rownames(scores), , drop = FALSE]

  # Una entrada por score, sin los sitios donde no se observó
  g <- NULL
  for (j in names(scores)) {
    g <- gstat::gstat(g       = g,
                      id      = j,
                      formula = stats::reformulate('1', response = j),
                      data    = scores[!is.na(scores[[j]]), j])
  }

  g
}

#' Translate fitted variogram models to the SpatFD representation
#'
#' Converts variogram models fitted outside the package into a common tabular
#' representation, the one [add_model()] registers in a `SpatFD` object.
#'
#' @param model Fitted models. Either a `gstat` object holding models (for
#'   example the output of [gstat::fit.lmc()]), or a named list of
#'   `variogramModel` objects (for example the outputs of
#'   [gstat::fit.variogram()]), one per score.
#' @param ... Currently unused.
#'
#' @return A data frame with one row per model structure and the columns
#' \describe{
#'   \item{`var`}{Id of the variogram the structure belongs to. For a `gstat`
#'     object fitted with cross-variograms this includes the cross ids.}
#'   \item{`model`}{Name of the model, as used in [models].}
#'   \item{`rango`}{Range parameter.}
#'   \item{`sigma`}{Partial sill.}
#'   \item{`alpha1`}{Additional shape parameter (`kappa` in \pkg{gstat}).}
#' }
#'
#' @seealso [togstat()], [add_model()], [availablemodels()]
#'
#' @examples
#' fitted <- list(
#'   "X.PC1" = gstat::vgm(psill = 2, model = "Exp", range = 300, nugget = 0.5),
#'   "X.PC2" = gstat::vgm(psill = 1, model = "Sph", range = 150)
#' )
#' toSpatFD(fitted)
#'
#' @export
toSpatFD <- function(model, ...) {
  UseMethod("toSpatFD")
}

#' @rdname toSpatFD
#' @export
toSpatFD.gstat <- function(model, ...) {
  gstat_names <- c(
    'Nug' = 'nugget',
    'Exp' = 'exponential',
    'Sph' = 'spherical',
    'Gau' = 'gaussian',
    'Exc' = 'stable',
    'Mat' = 'matern',
    'Ste' = "stein's matern",
    'Cir' = 'circular',
    'Lin' = 'linear',
    'Bes' = 'bessel',
    'Pen' = 'pentamodel',
    'Per' = 'periodic',
    'Wav' = 'wave',
    'Hol' = 'hole',
    'Log' = 'logarithmic',
    'Pow' = 'power',
    'Spl' = 'spline',
    'Leg' = 'legendre'
  )

  stopifnot('There are no models to extract from this object' = length(model$model) != 0)

  fitted <- model$model
  df <- do.call(rbind, fitted)
  df <- data.frame(var    = rep(names(fitted), vapply(fitted, nrow, integer(1))),
                   model  = unname(gstat_names[as.character(df$model)]),
                   rango  = df$range,
                   sigma  = df$psill,
                   alpha1 = df$kappa)
  df
}

#' @rdname toSpatFD
#' @export
toSpatFD.list <- function(model, ...){
  is_gstat <- vapply(model, inherits, logical(1), what = 'variogramModel')

  if (length(model) == 0 || !all(is_gstat) || is.null(names(model)) || any(names(model) == '')) {
    stop('The list must be a named list of gstat variogram models (variogramModel), ',
         'one per score', call. = FALSE)
  }

  toSpatFD(structure(list(model = model), class = c('gstat', 'list')))
}

#' @rdname toSpatFD
#' @export
toSpatFD.default <- function(model, ...) {
  stop("No SpatFD translator registered for objects of class '",
       class(model)[1], "'", call. = FALSE)
}
