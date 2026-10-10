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
#' site, when they have one.
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
#'     functional data object of \pkg{fda} with one curve per site.}
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
#' add_model(sfd, list(
#'   "sim.PC1" = gstat::vgm(psill = 30, model = "Exp", range = 0.5),
#'   "sim.PC2" = gstat::vgm(psill = 0.01, model = "Exp", range = 0.5)
#' ), name = "exponential")
#' add_model(sfd, list(
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
#' @importFrom stats predict
#' @export
leave_one_out <- function(sfd, design = NULL, variable = NULL){

  stopifnot('function only available for SpatFD objects' = inherits(sfd, 'SpatFD'))

  if (is.null(attr(sfd, 'args'))) {
    stop('The object does not keep the arguments of spatfd(). Create it again with spatfd()',
         call. = FALSE)
  }
  if (is.null(variable)) {
    variable <- var_names(sfd)[1]
    if (length(var_names(sfd)) > 1) message('Cross-validating the variable ', variable)
  }
  stopifnot(
    "'variable' must be the name of one variable of the object" =
      is.character(variable) && length(variable) == 1 && variable %in% var_names(sfd)
  )

  registry <- attr(sfd, 'models')
  ids      <- score_ids(sfd)
  own_ids  <- names(ids)[ids == variable]

  # =========================== COMBINACIONES ==================================
  # Por defecto: cada modelo que cubre la variable con cada método que admite
  # (lambda no admite modelos con variogramas cruzados)
  if (is.null(design)) {
    covers  <- vapply(as.list(registry), function(m) all(own_ids %in% names(m$variograms)), logical(1))
    crossed <- vapply(as.list(registry), function(m) !all(m$table$var %in% names(ids)), logical(1))
    order   <- rownames(sfd$models)

    design <- rbind(
      data.frame(model = order[covers[order]], method = rep('scores', sum(covers))),
      data.frame(model = order[covers[order] & !crossed[order]],
                 method = rep('lambda', sum(covers & !crossed)))
    )
  }
  stopifnot(
    "'design' must be a data frame with the columns 'model' and 'method'" =
      is.data.frame(design) && all(c('model', 'method') %in% colnames(design)),
    "'method' must be 'scores' or 'lambda'" = all(design$method %in% c('scores', 'lambda'))
  )
  if (nrow(design) == 0) {
    stop('There are no combinations to compare: no model covers the variable ', variable,
         '\nThe models of the object are:\n',
         paste(utils::capture.output(print(sfd$models)), collapse = '\n'), call. = FALSE)
  }
  unknown <- setdiff(design$model, ls(registry))
  if (length(unknown) > 0) {
    stop('Models not found in the object: ', paste(unknown, collapse = ', '),
         '\nThe models of the object are:\n',
         paste(utils::capture.output(print(sfd$models)), collapse = '\n'), call. = FALSE)
  }
  design <- unique(data.frame(model = as.character(design$model), method = as.character(design$method)))
  labels <- paste(design$model, design$method, sep = ' / ')

  # ============================ PREPARACIÓN ===================================
  observed <- attr(sfd, 'fd')[[variable]]                    # Curvas observadas (suavizadas)
  gram     <- fda::eval.penalty(observed$basis, fda::int2Lfd(0))   # Productos internos de la base
  sites    <- colnames(observed$coefs)
  coords   <- attr(sfd, 'data')$coords

  # Las demás variables no cambian de un sitio a otro: se reconstruyen una vez
  others <- lapply(setdiff(var_names(sfd), variable), function(v) refit(sfd, v))

  errors <- matrix(NA_real_, length(sites), nrow(design), dimnames = list(sites, labels))
  coefs  <- lapply(labels, function(l) observed$coefs * NA_real_)
  names(coefs) <- labels
  crossed <- stats::setNames(logical(nrow(design)), labels)

  # ========================= VALIDACIÓN CRUZADA ===============================
  for (site in sites) {

    # El análisis completo de la variable se repite sin el sitio
    fold <- Reduce(`+`, c(list(refit(sfd, variable, without = site)), others))
    for (model in unique(design$model)) add_model(fold, registry[[model]]$table, name = model)

    for (k in seq_len(nrow(design))) {
      # lambda solo predice una variable; scores usa todas las que cubre el modelo
      pred <- predict(fold, coords[site, , drop = FALSE], model = design$model[k],
                      method = design$method[k],
                      variable = if (design$method[k] == 'lambda') variable else NULL)
      if (!(variable %in% names(pred))) {
        stop("The model '", design$model[k], "' does not cover the variable ", variable,
             call. = FALSE)
      }

      coefs[[k]][, site] <- pred[[variable]]$pred$coefs
      crossed[k]         <- attr(pred, 'crossed')

      # Norma de la diferencia entre la curva observada y la predicha
      difference      <- observed$coefs[, site] - pred[[variable]]$pred$coefs[, 1]
      errors[site, k] <- sqrt(drop(crossprod(difference, gram %*% difference)))
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

  structure(list(summary = summary,
                 errors  = errors,
                 pred    = lapply(coefs, function(x) fda::fd(x, observed$basis, observed$fdnames))),
            class = 'SpatFDloo',
            variable = variable)
}

# FUNCIÓN DE IMPRESIÓN: Solo la tabla de comparación

#' @rdname leave_one_out
#' @param x An object of class `SpatFDloo`, as returned by `leave_one_out()`.
#' @param ... Ignored.
#' @export
print.SpatFDloo <- function(x, ...){
  cat('\n\n------ Leave-one-out cross-validation: ', attr(x, 'variable'), ' ------\n', sep = '')
  cat(nrow(x$errors), ' sites, error = distance between the observed and the predicted curve\n\n',
      sep = '')
  print(x$summary)
  cat('------------------------------------------------------\n\n')
  invisible(x)
}

# .............................. FUNCIONES INTERNAS ................................................

# REPETIR EL ANÁLISIS DE UNA VARIABLE: Vuelve a suavizar y a calcular los
# componentes principales de una variable con los mismos argumentos de
# spatfd(), opcionalmente sin uno de sus sitios. Los componentes quedan con el
# mismo signo que los del análisis original: el signo de un componente es
# arbitrario, pero los variogramas cruzados de un modelo dependen de él.
refit <- function(sfd, variable, without = NULL){
  args     <- attr(sfd, 'args')[[variable]]
  dat      <- attr(sfd, 'data')[[variable]]
  original <- attr(sfd, 'pca')[[variable]]
  scores   <- names(score_ids(sfd))[score_ids(sfd) == variable]

  if (!is.null(without)) dat <- dat[, colnames(dat) != without, drop = FALSE]

  new <- spatfd(dat, attr(sfd, 'data')$coords, t = args$t, basis = args$basis,
                nbasis = args$nbasis, lambda = args$lambda,
                nharm = length(scores), name = variable)

  # Signo de cada componente respecto al original, por su producto interno
  pca   <- attr(new, 'pca')[[variable]]
  gram  <- fda::eval.penalty(pca$harmonics$basis, fda::int2Lfd(0))
  signs <- sign(diag(crossprod(pca$harmonics$coefs, gram %*% original$harmonics$coefs)))
  signs[signs == 0] <- 1

  pca$harmonics$coefs <- sweep(pca$harmonics$coefs, 2, signs, `*`)
  pca$scores          <- sweep(pca$scores, 2, signs, `*`)
  attr(new, 'pca')[[variable]] <- pca

  new$scores[, scores] <- sweep(new$scores[, scores, drop = FALSE], 2, signs, `*`)
  new
}
