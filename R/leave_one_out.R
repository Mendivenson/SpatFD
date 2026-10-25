# =============================== leave_one_out() ==============================

#' Leave-one-out cross-validation of spatial functional models
#'
#' Compares models and prediction methods by leave-one-out cross-validation.
#' Each site of a variable is left out in turn, its curve is predicted from
#' the remaining sites with [predict()][predict.SpatFDmodel], and the
#' prediction is compared with the curve observed there.
#'
#' @param ... Objects of class `SpatFDmodel`, as returned by [spatfd_model()].
#'   They are labelled with the names of their arguments
#'   (`leave_one_out(a, other = b)`) or, without a name, with the expression
#'   given.
#' @param variable Character vector with the names of the variables to
#'   validate. Defaults to every variable of the models. A model is only
#'   validated on the variables it has.
#' @param method Prediction methods to compare: `"scores"`, `"lambda"` or
#'   both (the default). See [predict()][predict.SpatFDmodel].
#'
#' @details
#' **What is repeated for each site.** The whole analysis of the variable is
#' repeated without the site, with the same data, times and arguments given
#' to [spatfd()] and the same number of components of the model, so the curve
#' left out takes no part in the mean or in the components used to predict
#' it. The components obtained without the site are given the sign of the
#' ones of the whole analysis, so that the cross-variograms keep their
#' meaning. The variogram models are **not** fitted again: they are the ones
#' of the model. The other variables of the model keep their observation at
#' that site.
#'
#' **The error of a site.** The error of the prediction at a site is a curve
#' in time, \eqn{e(t) = \chi(t) - \hat\chi(t)}, the difference between the
#' observed (smoothed) curve and the predicted one: at each time, how far the
#' prediction was. It is summarised in a single number, the root mean square
#' error over time,
#' \deqn{RMSE = \sqrt{\frac{1}{|T|} \int_T e(t)^2 \, dt}}
#' computed as the mean over the times of the data (exactly the integral
#' when they are equally spaced). It is in the units of the data: a typical
#' distance between the predicted and the observed curve along the interval.
#' The predicted curves are kept in `pred`, so the error curves can be drawn.
#'
#' **Method `"lambda"` with several variables.** The curves of several
#' variables can only be combined over the times they share, so the method
#' `"lambda"` predicts, and its error is measured, only over the time
#' interval common to the variables of the model. When the variables span
#' different intervals, those errors are not directly comparable with the
#' ones of the method `"scores"`, which uses the whole interval of the
#' variable.
#'
#' @return An object of class `SpatFDloo`: a list with
#' \describe{
#'   \item{`summary`}{List with one table per variable, with one row per
#'     model and method: the `mean`, standard deviation (`sd`) and maximum
#'     (`max`) of the errors of the sites, from the smallest to the largest
#'     mean error.}
#'   \item{`errors`}{List with one table per variable, with the error of
#'     every site for every model and method.}
#'   \item{`pred`}{List with one element per variable and, inside, one per
#'     model and method (`<model>.<method>`): the curves predicted at each site
#'     when it was left out, with one row per time and one column per site.}
#' }
#'
#' @references
#' Bohórquez, M., Giraldo, R. and Mateu, J. (2017). Multivariate functional
#' random fields: prediction and optimal sampling. *Stochastic Environmental
#' Research and Risk Assessment*, 31, 53-70. \doi{10.1007/s00477-016-1266-y}
#'
#' @seealso [spatfd_model()], [predict()][predict.SpatFDmodel]
#'
#' @examples
#' set.seed(1)
#' dat <- sapply(1:10, function(i) rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1))
#' colnames(dat) <- paste0("s", 1:10)
#' coords <- cbind(x = runif(10), y = runif(10))
#' rownames(coords) <- colnames(dat)
#' sim <- spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "sim")
#'
#' exponential <- spatfd_model(sim, model = data.frame(var = c("sim.PC1", "sim.PC2"),
#'                                                     model = "exponential", rango = 0.5,
#'                                                     sigma = c(30, 0.1)))
#' spherical   <- spatfd_model(sim, model = data.frame(var = c("sim.PC1", "sim.PC2"),
#'                                                     model = "spherical", rango = 0.8,
#'                                                     sigma = c(30, 0.1)))
#'
#' loo <- leave_one_out(exponential, spherical)
#' loo
#' loo$summary$sim
#'
#' @export
leave_one_out <- function(...,                                # Modelos, ver spatfd_model()
                          variable = NULL,                    # Variables que se validan
                          method = c('scores', 'lambda')){    # Métodos de predicción

  models <- list(...)
  stopifnot(
    "At least one SpatFDmodel object must be given" = length(models) > 0,
    "All the objects in '...' must be of class SpatFDmodel" =
      all(vapply(models, inherits, logical(1), what = 'SpatFDmodel')),
    "'method' must be 'scores', 'lambda' or both" =
      is.character(method) && length(method) > 0 && all(method %in% c('scores', 'lambda'))
  )
  labels <- model_labels(models, as.list(substitute(list(...)))[-1])
  method <- unique(method)

  # VARIABLES: Todas las de los modelos, en orden de aparición
  available <- unique(unlist(lapply(models, function(m) names(get_model_data(m)))))
  if (is.null(variable)) variable <- available
  unknown <- setdiff(variable, available)
  if (length(unknown) > 0) {
    stop('Variables not in the models: ', paste(unknown, collapse = ', '), call. = F)
  }

  errors <- stats::setNames(vector('list', length(variable)), variable)
  pred   <- stats::setNames(lapply(variable, function(v) list()), variable)
  curves <- times <- pred

  for (k in seq_along(models)) {
    data  <- get_model_data(models[[k]])
    input <- get_compiled(models[[k]])$input

    for (v in intersect(variable, names(data))) {
      copy   <- data[[v]]
      coords <- get_coords(copy)

      for (site in rownames(coords)) {
        # MODELO SIN EL SITIO: La variable se rehace sin él y el modelo se arma
        # con los mismos variogramas
        fold_data      <- data
        fold_data[[v]] <- refit_without(copy, site)
        fold <- tryCatch(suppressMessages(do.call(spatfd_model, c(unname(fold_data), list(model = input)))),
                         error = function(e) {
                           stop("Without the site '", site, "' of '", v, "': ", conditionMessage(e), call. = F)
                         })

        for (m in method) {
          p <- suppressMessages(predict(fold, coords[site, , drop = F], variable = v, method = m))

          # Tiempos de la predicción: los de los datos, dentro del intervalo predicho
          range <- attr(p, 'internal')$pred[[v]]$basis$rangeval
          t     <- get_internal(fold_data[[v]], 't')
          t     <- t[t >= range[1] & t <= range[2]]

          # ERROR: La curva de error en el tiempo y su resumen (RMSE)
          e <- error_curve(get_fd(copy)[which(rownames(coords) == site)], p[[v]]$pred, t)
          errors[[v]] <- rbind(errors[[v]],
                               data.frame(site = site, model = labels[k], method = m, error = sqrt(mean(e^2))))

          # Curva predicha y curva de error en el sitio, como una columna más de
          # su combinación de modelo y método
          key <- paste(labels[k], m, sep = '.')
          pred[[v]][[key]]   <- cbind(pred[[v]][[key]], matrix(p[[v]]$pred[, 1], dimnames = list(rownames(p[[v]]$pred), site)))
          curves[[v]][[key]] <- cbind(curves[[v]][[key]], matrix(e, dimnames = list(NULL, site)))
          times[[v]][[key]]  <- t
        }
      }
    }
  }

  # Las curvas de error solo se usan para graficar, ver plot.SpatFDloo()
  structure(list(summary = lapply(errors, summary_table), errors = errors, pred = pred),
            class    = 'SpatFDloo',
            internal = list(curves = curves, times = times,
                            template = stats::setNames(lapply(variable, function(v) {
                              copies <- lapply(models, function(m) get_model_data(m)[[v]])
                              get_internal(Find(Negate(is.null), copies), 'time')
                            }), variable)))
}
