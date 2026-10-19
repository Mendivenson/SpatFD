# ================================ spatfd() ====================================

#' Create a SpatFD object: functional smoothing and functional PCA
#'
#' Wraps the core functional data workflow of the \pkg{fda} package
#' (basis smoothing followed by functional principal component analysis) into a
#' single object of class `SpatFD`. The scores of the principal components are
#' the quantities that are later modelled spatially.
#'
#' @param dat Numeric matrix of observations. Each column is a site (station)
#'   and each row is a measurement time. Column names identify the sites.
#' @param coords Numeric matrix of site coordinates. Row names must contain the
#'   column names of `dat`.
#' @param t Numeric vector with the times at which the rows of `dat` were
#'   measured. Must have unique values and length `nrow(dat)`.
#' @param basis Basis used for smoothing, either `"bsplines"` or `"fourier"`.
#' @param nbasis Number of basis functions.
#' @param lambda Non-negative smoothing parameter.
#' @param nharm Maximum number of principal components to compute. Cannot be
#'   greater than `nbasis`.
#' @param var Optional proportion of variance in `(0, 1]`. If given, the
#'   smallest number of principal components whose cumulative explained
#'   variance exceeds `var` is kept. If that proportion is not reached with
#'   `nharm` components, `nharm` components are kept and a message is shown.
#' @param name Name of the functional variable. Defaults to the expression
#'   passed as `dat`.
#'
#' @details
#' The functional representation of the data (smoothed curves and principal
#' components) is kept inside the object, hidden from the user: it is used by
#' [plot()][plot.SpatFD], [predict()][predict.SpatFD] and [leave_one_out()],
#' which always return plain matrices and data frames.
#'
#' @return An object of class `SpatFD`: a list with the elements
#' \describe{
#'   \item{`scores`}{Matrix with one row per site, containing the coordinates
#'     followed by the principal component scores. Score columns are named
#'     `<name>.PC1`, `<name>.PC2`, ...}
#'   \item{`variance`}{Data frame with one row per variable: basis type,
#'     number of principal components and proportion of explained variance.}
#'   \item{`models`}{Summary table of the variogram models registered in the
#'     object, one row per model. It is empty until a model is added with
#'     [add_model()].}
#' }
#'
#' @section Theory:
#' Each column of `dat` is a curve \eqn{\chi_{s_i}(t)} observed at a finite
#' set of points. It is first smoothed with `nbasis` basis functions
#' \eqn{B_l(t)}, by penalised least squares with smoothing parameter `lambda`:
#' \deqn{\chi_{s_i}(t) \approx \sum_{l=1}^{p} c_{il} B_l(t)}
#' The smoothed curves are then represented with their first \eqn{K}
#' functional principal components (Karhunen-Loève expansion):
#' \deqn{\chi_s(t) \approx \mu(t) + \sum_{k=1}^{K} f_k(s) \xi_k(t)}
#' The eigenfunctions \eqn{\xi_k(t)} are common to all sites, so the spatial
#' information is carried by the scores \eqn{f_k(s_i)}. These are the values
#' stored in `scores`, and the ones modelled spatially afterwards (Bohorquez,
#' Giraldo and Mateu, 2016). \eqn{K} is set by `nharm` or chosen from `var`.
#' See [SpatFD-package] for the complete picture.
#'
#' @references
#' Bohorquez, M., Giraldo, R. and Mateu, J. (2016). Optimal sampling for
#' spatial prediction of functional data. *Statistical Methods &
#' Applications*, 25(1), 39-54. \doi{10.1007/s10260-015-0340-9}
#'
#' Ramsay, J. O. and Silverman, B. W. (2005). *Functional Data Analysis*
#' (2nd ed.). Springer. \doi{10.1007/b98888}
#'
#' @seealso [var_names()], [togstat()], [`+.SpatFD`]
#'
#' @examples
#' set.seed(1)
#' t <- 1:50
#' dat <- sapply(1:10, function(i) rnorm(1, 5) * sin(2 * pi * t / 50) + rnorm(50, sd = 0.1))
#' colnames(dat) <- paste0("s", 1:10)
#' coords <- cbind(x = runif(10), y = runif(10))
#' rownames(coords) <- colnames(dat)
#'
#' sfd <- spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "sim")
#' sfd
#' summary(sfd)
#'
#' @export
spatfd <- function(dat,                      # Matriz de datos (nombres de columnas = estaciones)
                   coords,                   # Matriz de coordenadas (nombres de filas = estaciones)
                   t = 1:nrow(dat),          # Tiempos de medición, uno por fila de dat
                   basis = 'bsplines',       # Base para el suavizamiento
                   nbasis = 91,              # Número de bases para el suavizamiento
                   lambda = 0.0002,          # Parámetro de suavizamiento
                   nharm = 20,               # Cantidad máxima de PC
                   var = NULL,               # Varianza deseada
                   name = NULL){             # Nombre de la variable

  # ========================= VERIFICACIÓN DE DATOS ============================
  stopifnot(
    # Matrix
    "'dat' must be a non-empty matrix" = is.matrix(dat) && nrow(dat) > 0 && ncol(dat) > 0,
    "'coords' must be a non-empty matrix" = is.matrix(coords) && nrow(coords) > 0 && ncol(coords) > 0,

    # Numeric
    "'t' must be a numeric vector" = is.numeric(t),
    "'lambda' must be numeric and greater than or equal to zero" = is.numeric(lambda) && lambda >= 0,
    "'nharm' must be a positive numeric value" = is.numeric(nharm) && nharm > 0,
    "'nbasis' must be a positive numeric value" = is.numeric(nbasis) && nbasis > 0,
    "'var' must be numeric and between 0 and 1 (exclusive of 0, inclusive of 1) or null" = is.null(var) || (var > 0 && var <= 1),

    # Strings
    "'name' must be a character or null" = is.null(name) || is.character(name),

    # Consistency in given arguments
    "'nharm' cannot be greater than 'nbasis'" = nharm <= nbasis,
    "Length of 't' must equal the number of rows in 'dat'" = length(t) == nrow(dat),
    "Column names of 'dat' must match row names in 'coords'" = all(colnames(dat) %in% rownames(coords)),
    "'t' must contain unique values" = !any(duplicated(t))
  )

  if (is.unsorted(t)) {
    warning("'t' is not sorted. If dat is sorted consider sorting it.", call. = FALSE)
  }

  # Cuando no se da 'name', la variable se llama como el argumento 'dat'
  if (is.null(name)) name <- deparse(substitute(dat))

  # ========================== SUAVIZAMIENTO ===================================
  hr    <- range(t)
  oplfd <- fda::vec2Lfd(bwtvec = c(1, ncol(dat)), rangeval = hr)      # Penalización
  fbasis <- switch(basis,
                   'bsplines' = fda::create.bspline.basis(hr, nbasis),
                   'fourier'  = fda::create.fourier.basis(hr, nbasis),
                   stop("'basis' must be either 'bsplines' or 'fourier'", call. = F))

  fd <- fda::smooth.basis(argvals = t, dat, fda::fdPar(fdobj = fbasis, Lfdobj = oplfd, lambda))$fd

  # ================================ PCA =======================================
  pca <- fda::pca.fd(fd, nharm = nharm)

  # SELECCIÓN AUTOMÁTICA: Con 'var', el mínimo número de componentes que
  # alcanza la varianza deseada
  if (!is.null(var)){
    selected_pc <- which(cumsum(pca$varprop) > var)
    if (length(selected_pc) > 0){
      # Solo se repite el análisis si se necesitan menos componentes que los ya calculados
      if (selected_pc[1] < nharm) pca <- fda::pca.fd(fd, nharm = selected_pc[1])
      message('\t\033[32mVariance of ', var * 100, '% reached with ',
              selected_pc[1], ' principal components.\033[0m', appendLF = T)
    } else {
      message('\t\033[033mWARNING: Desired variance not reached. \n\t\t A variance of ', round(sum(pca$varprop) * 100,2),
              '% was reached with ', nharm, ' principal components.\033[0m', appendLF = T)
    }
  }

  # =========================== OBJETO RESULTANTE ==============================
  scores <- pca$scores
  dimnames(scores) <- list(colnames(dat), paste0(name, '.PC', seq_len(ncol(scores))))
  scores <- cbind(coords[colnames(dat), , drop = FALSE], scores)

  variance <- data.frame('basis' = pca$harmonics$basis$type,
                         'no. PC' = length(pca$varprop),
                         'variance' = sum(pca$varprop), check.names = F)
  rownames(variance) <- name

  # Lo funcional queda oculto, ver utils.R. Los argumentos del suavizamiento
  # permiten repetir el análisis sin algún sitio, ver leave_one_out().
  variables <- list()
  variables[[name]] <- list(data = dat, fd = fd, pca = pca,
                            args = list(t = t, basis = basis, nbasis = nbasis, lambda = lambda))

  new_spatfd(scores, variance, coords, variables)
}

# ======================== COMBINACIÓN DE VARIABLES ============================

#' Combine two SpatFD objects into a multivariate one
#'
#' Joins the variables of two `SpatFD` objects. Sites are merged as a full
#' join: a site observed for only one of the variables gets `NA` scores for
#' the other.
#'
#' @param x1,x2 Objects of class `SpatFD`, as returned by [spatfd()]. Their
#'   coordinates must have the same column names.
#'
#' @details
#' If a variable name is present in both objects, the one in `x2` is ignored
#' with a warning.
#'
#' The models registered in both objects are kept. If a model of `x2` has the
#' name of a model of `x1`, it is renamed sequentially (`"model 1"`,
#' `"model 2"`, ...) with a message.
#'
#' @section Theory:
#' With several functional variables, prediction is done by cokriging: the
#' scores of all the variables are modelled jointly, through their variograms
#' and cross-variograms, so that each variable borrows information from the
#' others. The variables do not need to be observed at the same sites
#' (Bohorquez, Giraldo and Mateu, 2017).
#'
#' @references
#' Bohorquez, M., Giraldo, R. and Mateu, J. (2017). Multivariate functional
#' random fields: prediction and optimal sampling. *Stochastic Environmental
#' Research and Risk Assessment*, 31, 53-70. \doi{10.1007/s00477-016-1266-y}
#'
#' @return An object of class `SpatFD` containing the variables of both
#'   objects. See [spatfd()] for its structure.
#'
#' @examples
#' set.seed(1)
#' t <- 1:50
#' coords <- cbind(x = runif(10), y = runif(10))
#' rownames(coords) <- paste0("s", 1:10)
#' sim <- function() {
#'   dat <- sapply(1:10, function(i) rnorm(1, 5) * sin(2 * pi * t / 50) + rnorm(50, sd = 0.1))
#'   colnames(dat) <- rownames(coords)
#'   dat
#' }
#'
#' a <- spatfd(sim(), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")
#' b <- spatfd(sim(), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "b")
#'
#' ab <- a + b
#' var_names(ab)
#'
#' @export
`+.SpatFD` <- function(x1, x2){

  if (!inherits(x2, 'SpatFD')){
    stop("Addition is only possible between two SpatFD objects.")
  }

  coords.x1 <- get_coords(x1)
  coords.x2 <- get_coords(x2)

  if (any(sort(colnames(coords.x1)) != sort(colnames(coords.x2)))){
    stop("Non compatible coordinates (Likely different names or dimensions)")
  }

  names_x2     <- var_names(x2)
  common_names <- intersect(var_names(x1), names_x2)

  if (length(common_names) > 0){
    warning('There are repeated variables between the two SpatFD objects. ',
            'Ignoring the repeated variables of the object ', deparse(substitute(x2)))
    names_x2 <- setdiff(names_x2, common_names)
  }

  # Coordenadas de todos los sitios
  new_places <- setdiff(rownames(coords.x2), rownames(coords.x1))
  coords     <- rbind(coords.x1, coords.x2[new_places, , drop = F])

  # Scores: outer join/full join por sitio. Se emparejan con match(): indexar
  # un data frame por nombre de fila hace coincidencias parciales (el sitio
  # 's1' tomaría los scores de 's10')
  only_scores <- function(x, keep) x$scores[, keep, drop = F]
  score.x1 <- only_scores(x1, !(colnames(x1$scores) %in% colnames(coords)))
  score.x2 <- only_scores(x2, names(score_ids(x2))[score_ids(x2) %in% names_x2])

  pts    <- sort(unique(c(rownames(score.x1), rownames(score.x2))))
  scores <- cbind(coords[pts, , drop = F],
                  score.x1[match(pts, rownames(score.x1)), , drop = F],
                  score.x2[match(pts, rownames(score.x2)), , drop = F])
  rownames(scores) <- pts

  variance  <- rbind(x1$variance, x2$variance[names_x2, , drop = F])
  variables <- c(attr(x1, 'internal')$variables, attr(x2, 'internal')$variables[names_x2])

  # Modelos registrados: primero los de x1 y luego los de x2, en su orden de registro
  models <- get_models(x1)
  for (model in names(get_models(x2))){
    new_name <- model
    if (model %in% names(models)){
      new_name <- next_model_name(names(models))
      message("Model '", model, "' of the second object renamed to '", new_name, "'")
    }
    models[[new_name]] <- get_models(x2)[[model]]
  }

  new_spatfd(scores, variance, coords, variables, models)
}
