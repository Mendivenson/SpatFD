# ================================ predict() ===================================

#' Spatial prediction of functional data
#'
#' Predicts the curves of one or several functional variables at unsampled
#' sites. The scores of the principal components are predicted by simple
#' kriging (or simple cokriging) and the curves are rebuilt from them. The
#' kriging system is solved inside the package, in matrix form.
#'
#' @param object An object of class `SpatFD`, as returned by [spatfd()].
#' @param newcoords Numeric matrix or data frame with the coordinates of the
#'   prediction sites, one row per site. Its columns must be the coordinates
#'   of `object`, matched by name when it has column names and by position
#'   otherwise. Row names, if any, identify the predicted curves.
#' @param model Name of one of the models registered in `object` with
#'   [add_model()]; `object$models` lists them. It can be omitted when
#'   the object has a single model. The model must have a variogram for every
#'   score of the variables in `variable`; cross-variograms between two scores
#'   are optional, see Details.
#' @param variable Character vector with the names of the variables to
#'   predict. Defaults to all the variables covered by the model. With
#'   several variables their scores are predicted jointly.
#' @param method Prediction method, a single one per call. Only `"scores"` is
#'   available for now.
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
#' `newcoords`.
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
#' @seealso [add_model()] to register models in the object.
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
#' plot(pred$sim$pred)
#' pred$sim$var
#'
#' @importFrom utils capture.output
#' @export
predict.SpatFD <- function(object, newcoords, model = NULL, variable = NULL, method = 'scores', ...){

  stopifnot("'method' must be a single string" = is.character(method) && length(method) == 1)
  method <- match.arg(method, c('scores', 'lambda'))
  if (method == 'lambda'){
    stop("The method 'lambda' is not available yet", call. = FALSE)
  }

  # Los scores predichos son un paso intermedio: al usuario solo se le entregan
  # las curvas y la varianza de predicción
  kriging <- krige_scores(object, newcoords, model, variable)
  sites   <- rownames(kriging$newcoords)

  # ====================== RECONSTRUCCIÓN DE LAS CURVAS ========================
  # Curva = media + combinación lineal de las funciones propias, todo sobre los
  # coeficientes de la base
  result <- lapply(names(kriging$pcs), function(v) {
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
  names(result) <- names(kriging$pcs)

  structure(result,
            class = 'SpatFDpred',
            newcoords = kriging$newcoords)
}

# FUNCIÓN DE IMPRESIÓN: Las curvas predichas (fd) no se imprimen, solo la
#                       varianza de predicción en cada sitio

#' @rdname predict.SpatFD
#' @param x An object of class `SpatFDpred`, as returned by `predict()`.
#' @export
print.SpatFDpred <- function(x, ...){
  cat('\n\n-------------------------------------------------\n')
  cat('Predicted curves at ', nrow(attr(x, 'newcoords')), ' sites for: ',
      paste(names(x), collapse = ', '), '\n', sep = '')
  cat('\nPrediction variance:\n')
  print(sapply(x, function(variable) variable$var))
  cat('-------------------------------------------------\n\n')
  invisible(x)
}

# .............................. FUNCIONES INTERNAS ................................................

# KRIGING DE LOS SCORES: Predice por kriging simple los scores de las variables
# pedidas en los sitios nuevos. Es el paso común a todo lo que necesite
# predicciones (predict(), validación cruzada, ...) y devuelve:
#   pred      -> Matriz de scores predichos (sitios nuevos x scores)
#   variance  -> Matriz con la varianza de predicción de cada score
#   pcs       -> Lista con los scores que pertenecen a cada variable
#   newcoords -> Coordenadas de los sitios nuevos, ya verificadas
krige_scores <- function(object, newcoords, model = NULL, variable = NULL){

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

  # ================================ KRIGING ===================================
  kriging <- simple_kriging(sites, values, score_covariances(compiled, ids), newcoords)
  dimnames(kriging$pred) <- dimnames(kriging$variance) <- list(rownames(newcoords), ids)

  list(pred = kriging$pred, variance = kriging$variance, pcs = pcs, newcoords = newcoords)
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
