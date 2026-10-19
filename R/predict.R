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
#' @param kriging Type of kriging: `"ordinary"` (the default), with an unknown
#'   constant mean, or `"simple"`, with a known mean of zero. See Details.
#' @param n Only when `newcoords` is a map: approximate number of points of
#'   the grid on which the curves are predicted.
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
#' By default the systems are solved by ordinary kriging: the mean is an
#' unknown constant, so the predictions do not depend on the curves being
#' centred, and models without a sill (`linear`, `power`, `logarithmic`,
#' `spline` and `legendre`) can be used. With `kriging = "simple"` the mean is
#' known and equal to zero, as the scores are centred; every model must then
#' have a sill.
#'
#' With the method `"lambda"` the curves of every variable in `variable` take
#' part in the prediction of each one (functional cokriging). The curves of two
#' variables are related through the cross-variograms between their scores; a
#' model without them predicts each variable on its own. The cross-variograms
#' between scores of the same variable play no role in this method.
#'
#' The curves of several variables can only be combined over the times they
#' share, so the `t` given to [spatfd()] must be on a common scale (for
#' example, hours since the same origin). When the variables span different
#' intervals, the method `"lambda"` predicts the curves over the common one
#' and says so with a message; the method `"scores"` is not affected.
#'
#' @return An object of class `SpatFDpred`: a list with one element per
#' predicted variable, named after it, each with
#' \describe{
#'   \item{`pred`}{Predicted curves, as a matrix with one row per time of the
#'     observed data (the `t` given to [spatfd()]) and one column per
#'     prediction site, the same layout as the data given to [spatfd()]. With
#'     the method `"lambda"` and variables observed over different intervals,
#'     only the times of the common interval.}
#'   \item{`var`}{Prediction variance of the curve at each prediction site,
#'     \eqn{E \| \chi_{s_0} - \hat{\chi}_{s_0} \|^2}.}
#' }
#' The list always has this form, also when a single variable is predicted.
#' The coordinates of the prediction sites are kept in the attribute
#' `newcoords`.
#'
#' The predictions can be plotted with [plot()][plot.SpatFDpred]: as curves
#' or, when the prediction was made over a map, as maps of the prediction and
#' of its variance.
#'
#' @section Theory:
#' Let \eqn{f} be the vector with the observed scores, of all the components
#' and variables involved, \eqn{\Sigma} its covariance matrix and \eqn{c_0} its
#' covariance with the score \eqn{f_k(s_0)} at the unsampled site \eqn{s_0}.
#' The scores are centred, so with simple kriging their mean is known and
#' equal to zero, and the predictor and its variance are
#' \deqn{\hat{f}_k(s_0) = \lambda^T f, \quad \Sigma \lambda = c_0, \qquad
#'   \sigma_k^2(s_0) = C_k(0) - c_0^T \lambda}
#' The covariances come from the fitted variograms, \eqn{C(h) = C(0) -
#' \gamma(h)}, where \eqn{C(0)} is the sum of the partial sills.
#'
#' Ordinary kriging takes the mean of each field as an unknown constant. The
#' weights of the field being predicted must add up to one and those of every
#' other field to zero (ordinary cokriging), which adds one Lagrange
#' multiplier \eqn{\mu} per field:
#' \deqn{\begin{pmatrix} \Sigma & F \\ F^T & 0 \end{pmatrix}
#'   \begin{pmatrix} \lambda \\ \mu \end{pmatrix} =
#'   \begin{pmatrix} c_0 \\ e \end{pmatrix}, \qquad
#'   \sigma_k^2(s_0) = C_k(0) - c_0^T \lambda - e^T \mu}
#' where \eqn{F} marks the observations of each field and \eqn{e} the field
#' being predicted. Under these constraints, adding a constant to a block of
#' covariances changes neither the weights nor the variance, so models without
#' a sill are used through \eqn{C(h) = -\gamma(h)} (Chilès and Delfiner, 2012).
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
#' of the observed curves, with one scalar weight per observed site and
#' variable (Giraldo, Delicado and Mateu, 2011; Bohorquez, Giraldo and Mateu,
#' 2017). For the variable \eqn{r},
#' \deqn{\hat{\chi}^r_{s_0}(t) = \mu^r(t) + \sum_{p=1}^{P} \sum_{i=1}^{n_p}
#'   \lambda^{rp}_i \left( \chi^p_{s_i}(t) - \mu^p(t) \right)}
#' The weights minimise \eqn{E \| \chi^r_{s_0} - \hat{\chi}^r_{s_0} \|^2} and
#' solve \eqn{\Omega \lambda = c}. \eqn{\Omega} has one block per pair of
#' variables, made of the covariances between curves
#' \deqn{E \langle \chi^p_{s}, \chi^q_{s'} \rangle = \sum_{k} \sum_{l}
#'   G^{pq}_{kl} \, C_{pk,ql}(\| s - s' \|), \qquad
#'   G^{pq}_{kl} = \int_T \xi^p_k(t) \, \xi^q_l(t) \, dt}
#' where \eqn{C_{pk,ql}} is the covariance between the scores \eqn{f^p_k} and
#' \eqn{f^q_l}, and \eqn{T} the interval common to the variables. For a
#' variable observed over all of \eqn{T}, \eqn{G^{pp}} is the identity, because
#' its eigenfunctions are orthonormal. \eqn{c} is built in the same way from
#' the covariances with \eqn{s_0}, and the prediction variance is
#' \eqn{E \langle \chi^r_{s_0}, \chi^r_{s_0} \rangle - c^T \lambda} (minus
#' \eqn{e^T \mu} with ordinary kriging), which uses the sills of the fitted
#' model. With ordinary kriging the weights of \eqn{r} add up to one and those
#' of every other variable to zero, so the means cancel out and the predictor
#' is the same whether the curves are centred or not:
#' \eqn{\hat{\chi}^r_{s_0}(t) = \sum_{p} \sum_{i} \lambda^{rp}_i
#' \chi^p_{s_i}(t)}. Unlike the method `"scores"`, the same weights are used
#' for all the scores and the curves are not truncated to the first \eqn{K}
#' components.
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
#' sfd <- add_model(sfd, list(
#'   "sim.PC1" = gstat::vgm(psill = 30, model = "Exp", range = 0.5),
#'   "sim.PC2" = gstat::vgm(psill = 0.01, model = "Exp", range = 0.5)
#' ), name = "exponential")
#' sfd$models
#'
#' pred <- predict(sfd, newcoords = cbind(x = c(0.2, 0.7), y = c(0.5, 0.5)), model = "exponential")
#' pred
#' pred$sim$var
#' plot(pred)
#'
#' # The same model with the method lambda, and with simple kriging
#' pred_lambda <- predict(sfd, newcoords = cbind(x = c(0.2, 0.7), y = c(0.5, 0.5)),
#'                        model = "exponential", method = "lambda")
#' plot(pred_lambda)
#' predict(sfd, newcoords = cbind(x = c(0.2, 0.7), y = c(0.5, 0.5)),
#'         model = "exponential", kriging = "simple")
#'
#' @export
predict.SpatFD <- function(object, newcoords, model = NULL, variable = NULL, method = 'scores',
                           kriging = 'ordinary', n = 5000, ...){

  stopifnot("'method' must be a single string" = is.character(method) && length(method) == 1,
            "'kriging' must be a single string" = is.character(kriging) && length(kriging) == 1)
  method  <- match.arg(method, c('scores', 'lambda'))
  kriging <- match.arg(kriging, c('ordinary', 'simple'))

  # MAPAS: Si en lugar de sitios se recibe un mapa, se predice sobre una grilla
  # regular de puntos dentro de él
  map <- NULL
  if (inherits(newcoords, c('sf', 'sfc'))) {
    map       <- newcoords
    newcoords <- map_grid(map, n)
  }

  # Cada método entrega las curvas predichas, evaluadas y como fd, y su
  # varianza; los pasos intermedios (scores predichos, pesos) no se le
  # entregan al usuario
  result <- switch(method,
                   'scores' = predict_scores(object, newcoords, model, variable, kriging),
                   'lambda' = predict_lambda(object, newcoords, model, variable, kriging))

  # El usuario recibe las curvas evaluadas en los tiempos observados. Las
  # curvas funcionales quedan ocultas, solo para graficar.
  variables <- names(result$curves)
  out <- lapply(result$curves, function(curve) list(pred = curve$pred, var = curve$var))

  structure(out,
            class     = 'SpatFDpred',
            newcoords = result$newcoords,
            method    = method,
            kriging   = kriging,
            crossed   = result$crossed,
            map       = map,
            internal  = list(pred     = lapply(result$curves, `[[`, 'fd'),
                             observed = lapply(stats::setNames(nm = variables), get_fd, sfd = object),
                             coords   = get_coords(object)))
}

# .............................. FUNCIONES INTERNAS ................................................

# MÉTODO SCORES: Kriging (o cokriging) de los scores y reconstrucción de las
# curvas como media + combinación lineal de las funciones propias, todo sobre
# los coeficientes de la base.
#
# Los dos métodos devuelven lo mismo, que es lo único que usa predict():
#   curves    -> Lista por variable con las curvas como fd (fd) y evaluadas en
#                los tiempos observados (pred), y su varianza (var)
#   newcoords -> Coordenadas de los sitios nuevos, ya verificadas
#   crossed   -> Si se usó algún variograma cruzado
predict_scores <- function(object, newcoords, model, variable, kriging = 'ordinary'){
  kriging <- krige_scores(object, newcoords, model, variable, kriging)
  sites   <- rownames(kriging$newcoords)

  curves <- lapply(names(kriging$pcs), function(v) {
    pca <- get_pca(object, v)
    pcs <- kriging$pcs[[v]]

    coefs <- pca$harmonics$coefs %*% t(kriging$pred[, pcs, drop = FALSE])
    coefs <- coefs + as.vector(pca$meanfd$coefs)
    colnames(coefs) <- sites

    fdnames      <- get_fd(object, v)$fdnames
    fdnames[[2]] <- sites

    fd <- fda::fd(coefs, pca$harmonics$basis, fdnames)

    list(fd   = fd,
         pred = evaluate_curves(fd, get_variable(object, v)$args$t),
         var  = rowSums(kriging$variance[, pcs, drop = FALSE]))
  })
  names(curves) <- names(kriging$pcs)

  list(curves = curves, newcoords = kriging$newcoords, crossed = kriging$crossed)
}

# MÉTODO LAMBDA (COKRIGING FUNCIONAL): La curva predicha de cada variable es
# su media más una combinación lineal de las curvas observadas, ya centradas,
# de todas las variables. Hay un peso por sitio observado y variable, igual
# para todos los tiempos.
#
# El sistema tiene un bloque por pareja de variables (p, q), con la covarianza
# entre sus curvas
#
#     Omega_pq(h) = suma_k suma_l  G_pq[k, l] * C_{pk, ql}(h)
#
# donde C_{pk, ql} es la covarianza entre los scores pk y ql, y G_pq el producto
# interno entre las funciones propias de las dos variables sobre el intervalo
# común de tiempo. Dentro de una variable observada en todo ese intervalo,
# G_pp es la identidad: solo cuenta la covarianza de cada score consigo mismo.
#
# Con kriging ordinario los pesos de la variable que se predice suman 1 y los de
# cada una de las demás suman 0. Entonces las medias se cancelan y la
# predicción es la misma con las curvas centradas o sin centrar.
predict_lambda <- function(object, newcoords, model, variable, kriging = 'ordinary'){
  setup     <- prediction_setup(object, newcoords, model, variable)
  variables <- names(setup$pcs)
  ids       <- unlist(setup$pcs, use.names = FALSE)
  ordinary  <- kriging == 'ordinary'
  covs      <- score_covariances(setup$compiled, ids, ordinary)
  newcoords <- setup$newcoords

  # ========================== INTERVALO COMÚN =================================
  # Las curvas de varias variables solo se pueden combinar en los tiempos que
  # comparten. Con una sola variable es su propio intervalo.
  ranges <- vapply(variables, function(v) get_fd(object, v)$basis$rangeval, numeric(2))
  domain <- c(max(ranges[1, ]), min(ranges[2, ]))
  own    <- ranges[1, ] == domain[1] & ranges[2, ] == domain[2]    # Variables observadas en todo el intervalo

  if (domain[1] >= domain[2]) {
    stop("The variables have no time in common, so the method 'lambda' cannot combine them. ",
         "The 't' given to spatfd() must be on a common scale", call. = FALSE)
  }
  if (!all(own)) {
    message("Method 'lambda': the curves are predicted over the common time interval [",
            format(domain[1]), ", ", format(domain[2]), "]")
  }

  # ===================== COVARIANZA ENTRE CURVAS ==============================
  inner <- function(p, q) {
    if (p == q && own[[p]]) return(diag(length(setup$pcs[[p]])))   # Funciones propias ortonormales
    inner_product(get_pca(object, p)$harmonics, get_pca(object, q)$harmonics, domain)
  }

  curve_covariance <- function(p, q) {
    G     <- inner(p, q)
    terms <- list()
    for (k in seq_along(setup$pcs[[p]])) {
      for (l in seq_along(setup$pcs[[q]])) {
        covariance <- covs[[setup$pcs[[p]][k], setup$pcs[[q]][l]]]
        if (!is.null(covariance) && G[k, l] != 0) terms <- c(terms, list(list(G = G[k, l], C = covariance)))
      }
    }
    if (length(terms) == 0) return(NULL)                # Variables sin covarianza cruzada
    function(h) Reduce(`+`, lapply(terms, function(term) term$G * term$C(h)))
  }

  P     <- length(variables)
  omega <- matrix(vector('list', P * P), P, P, dimnames = list(variables, variables))
  for (p in variables) for (q in variables) omega[p, q] <- list(curve_covariance(p, q))

  # ¿Se usó alguna covarianza cruzada entre variables?
  crossed <- !all(vapply(omega[row(omega) != col(omega)], is.null, logical(1)))

  # ============================ SISTEMA =======================================
  # Los scores de una variable comparten sitios: los de su primer score
  places <- lapply(variables, function(v) setup$sites[[match(setup$pcs[[v]][1], ids)]])
  names(places) <- variables
  rows <- split(seq_len(sum(vapply(places, nrow, integer(1)))),
                rep(variables, vapply(places, nrow, integer(1))))[variables]

  evaluate <- function(covariance, a, b) {
    d <- cross_dist(a, b)
    if (is.null(covariance)) return(matrix(0, nrow(d), ncol(d)))
    matrix(covariance(d), nrow(d), ncol(d))
  }

  Sigma <- matrix(0, length(unlist(rows)), length(unlist(rows)))
  for (p in variables) {
    for (q in variables) Sigma[rows[[p]], rows[[q]]] <- evaluate(omega[[p, q]], places[[p]], places[[q]])
  }
  if (!ordinary) {
    R <- tryCatch(chol(Sigma), error = function(e) {
      stop('The covariance matrix of the curves is not positive definite. ',
           'Check for duplicated sites or for an invalid model', call. = FALSE)
    })
  }

  # ============================ CURVAS ========================================
  # Curvas observadas centradas, evaluadas en los tiempos dados
  centred <- function(p, t) {
    observed <- fda::eval.fd(t, get_fd(object, p))
    colnames(observed) <- colnames(get_fd(object, p)$coefs)
    observed[, rownames(places[[p]]), drop = FALSE] - as.vector(fda::eval.fd(t, get_pca(object, p)$meanfd))
  }

  curves <- lapply(variables, function(r) {
    args <- get_variable(object, r)$args
    t    <- args$t[args$t >= domain[1] & args$t <= domain[2]]

    # Pesos de todos los sitios observados para cada sitio nuevo
    C0   <- do.call(rbind, lapply(variables, function(p) evaluate(omega[[p, r]], places[[p]], newcoords)))
    sill <- omega[[r, r]](0)

    if (ordinary) {
      # Una restricción por variable: los pesos de r suman 1, los de las demás 0
      targets  <- lapply(variables, function(p) if (p == r) seq_len(nrow(newcoords)) else integer(0))
      solved   <- constrained_solve(Sigma, C0, rows, targets, sill)
      weights  <- solved$weights
      variance <- solved$variance
    } else {
      weights  <- backsolve(R, backsolve(R, C0, transpose = TRUE))
      variance <- sill - colSums(C0 * weights)
    }
    names(variance) <- rownames(newcoords)

    pred <- as.vector(fda::eval.fd(t, get_pca(object, r)$meanfd)) +
      Reduce(`+`, lapply(variables, function(p) centred(p, t) %*% weights[rows[[p]], , drop = FALSE]))
    dimnames(pred) <- list(t, rownames(newcoords))

    # Versión funcional, solo para graficar: en la base de la variable si se
    # predijo en todo su intervalo, en B-splines sobre el intervalo común si no
    basis <- if (own[[r]]) get_fd(object, r)$basis else fda::create.bspline.basis(domain, args$nbasis)
    fd    <- fda::smooth.basis(t, pred, basis)$fd
    fd$fdnames[[2]]  <- rownames(newcoords)
    colnames(fd$coefs) <- rownames(newcoords)

    list(fd = fd, pred = pred, var = variance)
  })
  names(curves) <- variables

  list(curves = curves, newcoords = newcoords, crossed = crossed)
}

# PREPARACIÓN COMÚN A LOS DOS MÉTODOS: Verifica los argumentos y deja listo lo
# que se necesita para predecir:
#   model, compiled -> Nombre del modelo escogido y su versión compilada
#   pcs             -> Lista con los scores que pertenecen a cada variable
#   sites, values   -> Sitios y valores observados de cada score
#   newcoords       -> Coordenadas de los sitios nuevos, ya verificadas
prediction_setup <- function(object, newcoords, model = NULL, variable = NULL){

  # ============================ MODELO A USAR =================================
  # Los modelos ya están compilados dentro del objeto, ver add_model(). En
  # caso de error se muestra la tabla resumen de los modelos registrados.
  available <- names(get_models(object))

  if (length(available) == 0) {
    stop('The object has no models. Register one with add_model()', call. = FALSE)
  }
  if (is.null(model)) {
    if (length(available) > 1) {
      stop("The object has several models, choose one with 'model':\n", models_message(object),
           call. = FALSE)
    }
    model <- available
  }
  if (!is.character(model) || length(model) != 1 || !(model %in% available)) {
    stop("'model' must be the name of one of the models of the object:\n", models_message(object),
         call. = FALSE)
  }
  compiled <- get_models(object)[[model]]

  # Scores del objeto (un campo = un score) y variable a la que pertenecen
  all_ids <- score_ids(object)

  if (is.null(variable)) {
    variable <- unique(all_ids[intersect(names(all_ids), names(compiled$variograms))])
  }

  coords <- get_coords(object)

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
         '\nThe models of the object are:\n', models_message(object), call. = FALSE)
  }

  # Sitios con score faltante se descartan campo a campo: las variables pueden
  # estar observadas en conjuntos de sitios distintos
  scores   <- object$scores
  observed <- lapply(ids, function(id) !is.na(scores[, id]))
  sites    <- lapply(observed, function(obs) scores[obs, colnames(coords), drop = FALSE])
  values   <- Map(function(id, obs) scores[obs, id], ids, observed)

  list(model = model, compiled = compiled, pcs = pcs, sites = sites, values = values,
       newcoords = newcoords)
}

# KRIGING DE LOS SCORES: Predice por kriging ordinario o simple los scores de
# las variables pedidas en los sitios nuevos. Es el paso común a todo lo que
# necesite los scores predichos (método scores, validación cruzada, ...) y
# devuelve:
#   pred      -> Matriz de scores predichos (sitios nuevos x scores)
#   variance  -> Matriz con la varianza de predicción de cada score
#   pcs       -> Lista con los scores que pertenecen a cada variable
#   newcoords -> Coordenadas de los sitios nuevos, ya verificadas
#   crossed   -> Si se usó algún variograma cruzado entre los scores predichos
krige_scores <- function(object, newcoords, model = NULL, variable = NULL, kriging = 'ordinary'){
  setup    <- prediction_setup(object, newcoords, model, variable)
  ids      <- unlist(setup$pcs, use.names = FALSE)
  ordinary <- kriging == 'ordinary'

  covs    <- score_covariances(setup$compiled, ids, ordinary)
  kriging <- solve_kriging(setup$sites, setup$values, covs, setup$newcoords, ordinary)
  dimnames(kriging$pred) <- dimnames(kriging$variance) <- list(rownames(setup$newcoords), ids)

  # ¿Se usó algún variograma cruzado entre los scores predichos?
  crossed <- !all(vapply(covs[upper.tri(covs)], is.null, logical(1)))

  list(pred = kriging$pred, variance = kriging$variance, pcs = setup$pcs,
       newcoords = setup$newcoords, crossed = crossed)
}
