#' Translate fitted variogram models to the SpatFD representation
#'
#' Converts variogram models fitted outside the package into a common tabular
#' representation, the one [add_model()] registers in a `SpatFD` object.
#'
#' @param model Fitted models. Either a `gstat` object holding models (for
#'   example the output of [gstat::fit.lmc()]), or a named list of
#'   `variogramModel` objects (for example the outputs of
#'   [gstat::fit.variogram()]), one per score.
#' @param sfd Currently unused.
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

# MÉTODO IMPLEMENTADO PARA LOS MODELOS PROVENIENTES DE GSTAT.

#' @rdname toSpatFD
#' @export
toSpatFD.gstat <- function(model, sfd, ...) {
  model.trans <- c(
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

  stopifnot(
    'There are no models to extract from this object' = length(model$model) != 0
  )

  model <- model$model
  vars  <- rep(names(model), times = unlist(lapply(model, nrow)))
  model <- model |>
    do.call(what = rbind)
  model <- cbind('var' = rownames(model), model)
  model <- model[,c('var', 'model', 'range', 'psill', 'kappa')]
  colnames(model) <- c('var', 'model', 'rango', 'sigma', 'alpha1')
  rownames(model) <- NULL

  model[,"model"] = model.trans[model$model]
  model$var = vars

  return(model)
}

# CUANDO SE TENGAN LISTAS DE MODELOS EN EL CASO UNIVARIADO

#' @rdname toSpatFD
#' @export
toSpatFD.list <- function(model,...){

  # Modelos GSTAT
  if (all(c('variogramModel', 'data.frame') == lapply(model, class) |>
      unlist() |>
      unique())){
        x = list()
        x$model = model
        class(x) = c('gstat', class(x))
  }

  toSpatFD(x)
}

#' @rdname toSpatFD
#' @export
toSpatFD.default <- function(model, ...) {
  stop("No SpatFD translator registered for objects of class '",
       class(model)[1], "'", call. = FALSE)
}
