# ============================== leave_one_out() ===============================

#' Leave-one-out cross-validation of models and prediction methods
#'
#' Compares variogram models and prediction methods by leave-one-out
#' cross-validation. Each site of a variable is left out in turn, its curve is
#' predicted from the remaining sites with [predict()][predict.SpatFD], and
#' the prediction is compared with the curve that was observed there.
#'
#' @param sfd An object of class `SpatFD` with at least one model registered,
#'   see [add_model()].
#' @param design Data frame with the combinations to compare, one per row,
#'   with the columns `model` (name of a model registered in `sfd`) and
#'   `method` (`"scores"` or `"lambda"`). By default, every model that covers
#'   the variable is combined with every method it can be used with.
#' @param variable Name of the variable whose curves are predicted. Defaults
#'   to the first variable of the object, see [var_names()].
#' @param kriging Type of kriging used in every combination: `"ordinary"` (the
#'   default) or `"simple"`, see [predict()][predict.SpatFD]. With
#'   `"simple"`, the default design leaves out the models without a sill.
#'
#' @details
#' For every site left out, the whole analysis of the variable is repeated
#' without it, with the same arguments given to [spatfd()] and the same number
#' of principal components: the curve left out takes no part in the mean
#' function or in the principal components used to predict it. The variogram
#' models are **not** fitted again; they are the ones registered in `sfd`.
#'
#' The principal components obtained without a site are given the same sign
#' as the ones of the full analysis, so that the cross-variograms of a model
#' keep their meaning.
#'
#' With a model that covers several variables (cokriging), only the curve of
#' `variable` is left out: the other variables keep their observation at that
#' site, when they have one. Both methods use every variable of the model.
#'
#' With the method `"lambda"` and variables observed over different time
#' intervals, the curves are predicted, and the errors measured, only over the
#' common interval (see [predict()][predict.SpatFD]). Those errors are then
#' not directly comparable with the ones of the method `"scores"`.
#'
#' @return An object of class `SpatFDloo`: a list with
#' \describe{
#'   \item{`summary`}{Data frame with one row per combination, from the
#'     smallest to the largest mean error: `model`, `method`, `crossed`
#'     (whether cross-variograms were used), `mean` and `sd` (mean and
#'     standard deviation of the errors over the sites) and `max` (largest
#'     error).}
#'   \item{`errors`}{Matrix with the error of every site (rows) for every
#'     combination (columns).}
#'   \item{`pred`}{List with the curves predicted for every combination, as a
#'     matrix with one row per time of the observed data and one column per
#'     site, the same layout as the data given to [spatfd()].}
#' }
#'
#' @section Theory:
#' The error at the site \eqn{s_i} is the distance between the observed
#' (smoothed) curve and the curve predicted without it,
#' \deqn{\| \chi_{s_i} - \hat{\chi}_{s_i}^{(-i)} \| = \left( \int_T \left(
#'   \chi_{s_i}(t) - \hat{\chi}_{s_i}^{(-i)}(t) \right)^2 dt \right)^{1/2}}
#' A combination of model and method predicts better the smaller its mean
#' error is (Bohorquez, Giraldo and Mateu, 2017).
#'
#' @references
#' Bohorquez, M., Giraldo, R. and Mateu, J. (2017). Multivariate functional
#' random fields: prediction and optimal sampling. *Stochastic Environmental
#' Research and Risk Assessment*, 31, 53-70. \doi{10.1007/s00477-016-1266-y}
#'
#' @seealso [add_model()], [predict()][predict.SpatFD]
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
#' ), name = "exponential")
#' sfd <- add_model(sfd, list(
#'   "sim.PC1" = gstat::vgm(psill = 30, model = "Sph", range = 0.8),
#'   "sim.PC2" = gstat::vgm(psill = 0.01, model = "Sph", range = 0.8)
#' ), name = "spherical")
#'
#' # Every model with every method
#' loo <- leave_one_out(sfd)
#' loo
#'
#' # Or only some combinations
#' leave_one_out(sfd, design = data.frame(model  = c("exponential", "spherical"),
#'                                        method = c("scores", "scores")))
#'
#' @export
leave_one_out <- function(sfd, design = NULL, variable = NULL, kriging = 'ordinary'){

  stopifnot('function only available for SpatFD objects' = inherits(sfd, 'SpatFD'),
            "'kriging' must be a single string" = is.character(kriging) && length(kriging) == 1)
  kriging <- match.arg(kriging, c('ordinary', 'simple'))

  if (is.null(variable)) {
    variable <- var_names(sfd)[1]
    if (length(var_names(sfd)) > 1) message('Cross-validating the variable ', variable)
  }
  stopifnot(
    "'variable' must be the name of one variable of the object" =
      is.character(variable) && length(variable) == 1 && variable %in% var_names(sfd)
  )

  registry <- get_models(sfd)
  ids      <- score_ids(sfd)
  own_ids  <- names(ids)[ids == variable]

  # =========================== COMBINACIONES ==================================
  # Por defecto: cada modelo que cubre la variable con cada método. Con kriging
  # simple quedan por fuera los modelos sin meseta.
  if (is.null(design)) {
    usable <- vapply(registry, function(m) {
      all(own_ids %in% names(m$variograms)) &&
        (kriging == 'ordinary' || !any(m$table$model %in% unbounded_models))
    }, logical(1))
    covers <- as.character(names(registry))[usable]
    design <- data.frame(model  = rep(covers, 2),
                         method = rep(c('scores', 'lambda'), each = length(covers)))
  }
  stopifnot(
    "'design' must be a data frame with the columns 'model' and 'method'" =
      is.data.frame(design) && all(c('model', 'method') %in% colnames(design)),
    "'method' must be 'scores' or 'lambda'" = all(design$method %in% c('scores', 'lambda'))
  )
  if (nrow(design) == 0) {
    stop('There are no combinations to compare: no model covers the variable ', variable,
         '\nThe models of the object are:\n', models_message(sfd), call. = FALSE)
  }
  unknown <- setdiff(design$model, names(registry))
  if (length(unknown) > 0) {
    stop('Models not found in the object: ', paste(unknown, collapse = ', '),
         '\nThe models of the object are:\n', models_message(sfd), call. = FALSE)
  }
  design <- unique(data.frame(model = as.character(design$model), method = as.character(design$method)))
  labels <- paste(design$model, design$method, sep = ' / ')

  # ============================ PREPARACIÓN ===================================
  observed <- get_fd(sfd, variable)                               # Curvas observadas (suavizadas)
  gram     <- fda::eval.penalty(observed$basis, fda::int2Lfd(0))   # Productos internos de la base
  sites    <- colnames(observed$coefs)
  times    <- get_variable(sfd, variable)$args$t
  coords   <- get_coords(sfd)

  # Las demás variables no cambian de un sitio a otro: se reconstruyen una vez
  others <- lapply(setdiff(var_names(sfd), variable), function(v) refit(sfd, v))

  errors <- matrix(NA_real_, length(sites), nrow(design), dimnames = list(sites, labels))
  curves <- lapply(labels, function(l) matrix(NA_real_, length(times), length(sites),
                                              dimnames = list(times, sites)))
  names(curves) <- labels
  crossed <- stats::setNames(logical(nrow(design)), labels)

  # ========================= VALIDACIÓN CRUZADA ===============================
  for (site in sites) {

    # El análisis completo de la variable se repite sin el sitio
    fold <- Reduce(`+`, c(list(refit(sfd, variable, without = site)), others))
    for (model in unique(design$model)) fold <- add_model(fold, registry[[model]]$table, name = model)

    for (k in seq_len(nrow(design))) {
      # Los dos métodos usan todas las variables que cubre el modelo
      pred <- suppressMessages(predict(fold, coords[site, , drop = FALSE], model = design$model[k],
                                       method = design$method[k], kriging = kriging))
      if (!(variable %in% names(pred))) {
        stop("The model '", design$model[k], "' does not cover the variable ", variable,
             call. = FALSE)
      }

      # Con lambda y variables en intervalos distintos solo se predice el común
      curve <- pred[[variable]]$pred
      curves[[k]][rownames(curve), site] <- curve[, 1]
      crossed[k] <- attr(pred, 'crossed')

      errors[site, k] <- curve_distance(observed[which(sites == site)], predicted_fd(pred, variable)[1], gram)
    }
  }

  # ============================== RESULTADOS ==================================
  summary <- data.frame(model   = design$model,
                        method  = design$method,
                        crossed = unname(crossed),
                        mean    = colMeans(errors),
                        sd      = apply(errors, 2, stats::sd),
                        max     = apply(errors, 2, max),
                        row.names = NULL)
  summary <- summary[order(summary$mean), ]
  rownames(summary) <- NULL

  structure(list(summary = summary, errors = errors, pred = curves),
            class = 'SpatFDloo',
            variable = variable,
            kriging = kriging)
}

# FUNCIÓN DE IMPRESIÓN: Solo la tabla de comparación

#' @rdname leave_one_out
#' @param x An object of class `SpatFDloo`, as returned by `leave_one_out()`.
#' @param ... Ignored.
#' @export
print.SpatFDloo <- function(x, ...){
  cat('\n\n------ Leave-one-out cross-validation: ', attr(x, 'variable'), ' ------\n', sep = '')
  cat(nrow(x$errors), ' sites, ', attr(x, 'kriging'), ' kriging, ',
      'error = distance between the observed and the predicted curve\n\n', sep = '')
  print(x$summary)
  cat('------------------------------------------------------\n\n')
  invisible(x)
}

# .............................. FUNCIONES INTERNAS ................................................

# DISTANCIA ENTRE DOS CURVAS: Norma L2 de su diferencia, sobre el intervalo de
# la curva predicha. Si están en la misma base basta con los coeficientes y la
# matriz de productos internos de la base (gram); si no, se integra.
curve_distance <- function(observed, predicted, gram){
  if (identical(observed$basis, predicted$basis)) {
    difference <- observed$coefs[, 1] - predicted$coefs[, 1]
    return(sqrt(drop(crossprod(difference, gram %*% difference))))
  }

  rng <- predicted$basis$rangeval
  squared <- inner_product(observed, observed, rng) - 2 * inner_product(observed, predicted, rng) +
    inner_product(predicted, predicted, rng)
  sqrt(max(0, drop(squared)))
}

# REPETIR EL ANÁLISIS DE UNA VARIABLE: Vuelve a suavizar y a calcular los
# componentes principales de una variable con los mismos argumentos de
# spatfd(), opcionalmente sin uno de sus sitios. Los componentes quedan con el
# mismo signo que los del análisis original: el signo de un componente es
# arbitrario, pero los variogramas cruzados de un modelo dependen de él.
refit <- function(sfd, variable, without = NULL){
  args     <- get_variable(sfd, variable)$args
  dat      <- get_variable(sfd, variable)$data
  original <- get_pca(sfd, variable)
  scores   <- names(score_ids(sfd))[score_ids(sfd) == variable]

  if (!is.null(without)) dat <- dat[, colnames(dat) != without, drop = FALSE]

  new <- spatfd(dat, get_coords(sfd), t = args$t, basis = args$basis,
                nbasis = args$nbasis, lambda = args$lambda,
                nharm = length(scores), name = variable)

  # Signo de cada componente respecto al original, por su producto interno
  pca   <- get_pca(new, variable)
  gram  <- fda::eval.penalty(pca$harmonics$basis, fda::int2Lfd(0))
  signs <- sign(diag(crossprod(pca$harmonics$coefs, gram %*% original$harmonics$coefs)))
  signs[signs == 0] <- 1

  pca$harmonics$coefs <- sweep(pca$harmonics$coefs, 2, signs, `*`)
  pca$scores          <- sweep(pca$scores, 2, signs, `*`)
  new <- set_pca(new, variable, pca)

  new$scores[, scores] <- sweep(new$scores[, scores, drop = FALSE], 2, signs, `*`)
  new
}
