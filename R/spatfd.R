# ================================ spatfd() ====================================

#' Create a SpatFD object: functional smoothing and functional PCA
#'
#' Wraps the core functional data workflow of the \pkg{fda} package (basis
#' smoothing followed by functional principal component analysis) into a
#' single object of class `SpatFD`, for one functional variable. The scores of
#' the principal components are the quantities that are later modelled
#' spatially.
#'
#' @param dat Data frame or numeric matrix of observations. Each row is a
#'   measurement time and each column a site (station), except for the column
#'   `time`. Column names identify the sites. Missing values are allowed, see
#'   Details.
#' @param coords Data frame or numeric matrix with the coordinates of the
#'   sites, one row per site. The sites are identified by the column `id` or,
#'   when `id` is `NULL`, by the row names.
#' @param time Name of the column of `dat` with the time of each row: numbers,
#'   dates (`Date`) or dates and times (`POSIXct`). If `dat` has no such
#'   column, the rows are taken as equally spaced times `1, ..., n` in their
#'   current order, and a message says so.
#' @param id Name of the column of `coords` that identifies the sites, or
#'   `NULL` to use the row names of `coords`. Two sites cannot have the same
#'   coordinates.
#' @param basis Basis used for smoothing, either `"bsplines"` or `"fourier"`.
#' @param nbasis Number of basis functions.
#' @param lambda Non-negative smoothing parameter.
#' @param nharm Maximum number of principal components to compute. Cannot be
#'   greater than `nbasis`.
#' @param var Optional proportion of variance in `(0, 1]`. If given, the
#'   smallest number of principal components whose cumulative explained
#'   variance exceeds `var` is kept. If that proportion is not reached with
#'   `nharm` components, `nharm` components are kept and a message is shown.
#' @param name Name of the functional variable. Defaults to the name of the
#'   object passed as `dat` (or the expression, if it is not a single
#'   object). It is kept in the attribute `name` and used to name the scores.
#'
#' @details
#' **Times.** Each row of `dat` is taken as measured at the time in its
#' `time` column, so the rows can come in any order: they are sorted by time.
#' For the smoothing, dates are converted to days since 1970-01-01 UTC, with
#' the time of day as a fraction of the day: a date is 00:00 UTC of that day.
#' Using UTC keeps the times increasing, with no hour repeated or skipped by
#' daylight saving time, and puts a variable with dates and one with dates
#' and times on the same axis. The times kept in the object are the ones
#' given, in their own time zone. Numeric times are used as they are, so different
#' variables must use the same scale. The unit matters for the smoothing: the
#' same `lambda` smooths differently in days than in hours.
#'
#' **Missing values.** Every curve is smoothed with all its available
#' observations: a site with missing values is smoothed on its own, over the
#' times where it was observed. Sites with no observations at all are
#' removed, with a message, and a warning lists the sites with more than 10%
#' of missing values. The 10% threshold is a heuristic: how well a curve is
#' recovered depends on how many observations are left compared with
#' `nbasis`, and on how they are spread (a long gap is worse than scattered
#' missing values). Rows with no observation at any site are dropped. The
#' percentage of missing values of each site is kept in the element
#' `missing`.
#'
#' The functional representation of the data (smoothed curves and principal
#' components) is kept inside the object, hidden from the user: it is used by
#' [plot()][plot.SpatFD].
#'
#' @return An object of class `SpatFD`: a list with the elements
#' \describe{
#'   \item{`scores`}{Matrix with one row per site, containing the coordinates
#'     followed by the principal component scores. Score columns are named
#'     `<name>.PC1`, `<name>.PC2`, ...}
#'   \item{`variance`}{Data frame with one row, named after the variable:
#'     basis type, number of principal components and proportion of
#'     explained variance.}
#'   \item{`missing`}{Named vector with the percentage of missing values of
#'     each site, over the rows kept.}
#' }
#' The name of the variable is kept in the attribute `name`:
#' `attr(x, "name")`.
#'
#' @section Theory:
#' Each column of `dat` is a curve \eqn{\chi_{s_i}(t)} observed at a finite
#' set of points. It is first smoothed with `nbasis` basis functions
#' \eqn{B_l(t)}:
#' \deqn{\chi_{s_i}(t) \approx \sum_{l=1}^{p} c_{il} B_l(t)}
#' The coefficients minimise the penalised least squares criterion
#' \deqn{\sum_{j} \left( y_{ij} - \chi_{s_i}(t_j) \right)^2 +
#'   \lambda \int_T \left( \chi_{s_i}''(t) \right)^2 dt}
#' over the times \eqn{t_j} where the site was observed, so the curves do not
#' need to share their times. The smoothed curves are then represented with
#' their first \eqn{K} functional principal components (Karhunen-Loève
#' expansion):
#' \deqn{\chi_s(t) \approx \mu(t) + \sum_{k=1}^{K} f_k(s) \xi_k(t)}
#' The eigenfunctions \eqn{\xi_k(t)} are common to all sites, so the spatial
#' information is carried by the scores \eqn{f_k(s_i)}. These are the values
#' stored in `scores`. \eqn{K} is set by `nharm` or chosen from `var`.
#'
#' @references
#' Bohórquez, M. P. (2025). *Geoestadística para datos espaciales,
#' espacio-temporales y funcionales* (1st ed.). Bogotá: Editorial UNAL.
#' ISBN 978-628-503-048-2.
#'
#' Bohórquez, M., Giraldo, R. and Mateu, J. (2016). Optimal sampling for
#' spatial prediction of functional data. *Statistical Methods &
#' Applications*, 25(1), 39-54. \doi{10.1007/s10260-015-0340-9}
#'
#' Bohórquez, M., Giraldo, R. and Mateu, J. (2017). Multivariate functional
#' random fields: prediction and optimal sampling. *Stochastic Environmental
#' Research and Risk Assessment*, 31, 53-70. \doi{10.1007/s00477-016-1266-y}
#'
#' Horváth, L. and Kokoszka, P. (2012). *Inference for Functional Data with
#' Applications*. Springer Series in Statistics. Springer.
#' \doi{10.1007/978-1-4614-3655-3}
#'
#' Ramsay, J. (2005). Functional Data Analysis. In B. S. Everitt and D. C.
#' Howell (eds.), *Encyclopedia of Statistics in Behavioral Science*.
#' \doi{10.1002/0470013192.bsa239}
#'
#' @seealso [plot()][plot.SpatFD], [summary()][SpatFD-print]
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(time = as.Date("2026-01-01") + 0:49)
#' for (s in paste0("s", 1:10)) {
#'   dat[[s]] <- rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1)
#' }
#' coords <- data.frame(station = paste0("s", 1:10), x = runif(10), y = runif(10))
#'
#' sfd <- spatfd(dat, coords, id = "station", basis = "fourier", nbasis = 11, nharm = 2, name = "sim")
#' sfd
#' summary(sfd)
#'
#' @export
spatfd <- function(dat,                      # Data frame/matriz de datos (nombres de columnas = estaciones)
                   coords,                   # Data frame/matriz de coordenadas
                   time = 'time',            # Columna de 'dat' con los tiempos
                   id = NULL,                # Columna de 'coords' con las estaciones (NULL = nombres de filas)
                   basis = 'bsplines',       # Base para el suavizamiento
                   nbasis = 91,              # Número de bases para el suavizamiento
                   lambda = 0.0002,          # Parámetro de suavizamiento
                   nharm = 20,               # Cantidad máxima de PC
                   var = NULL,               # Varianza deseada
                   name = NULL){             # Nombre de la variable

  # Cuando no se da 'name', la variable se llama como el argumento 'dat'. Va
  # antes de cualquier cambio a 'dat': después substitute() ya no da su nombre
  if (is.null(name)) name <- deparse(substitute(dat))

  # ========================= VERIFICACIÓN DE DATOS ============================
  stopifnot(
    # Data frame o matriz
    "'dat' must be a non-empty data frame or matrix" = (is.data.frame(dat) || is.matrix(dat)) && nrow(dat) > 0 && ncol(dat) > 0,
    "'coords' must be a non-empty data frame or matrix" = (is.data.frame(coords) || is.matrix(coords)) && nrow(coords) > 0 && ncol(coords) > 0,

    # Numeric
    "'lambda' must be numeric and greater than or equal to zero" = is.numeric(lambda) && lambda >= 0,
    "'nharm' must be a positive numeric value" = is.numeric(nharm) && nharm > 0,
    "'nbasis' must be a positive numeric value" = is.numeric(nbasis) && nbasis > 0,
    "'var' must be numeric and between 0 and 1 (exclusive of 0, inclusive of 1) or null" = is.null(var) || (var > 0 && var <= 1),

    # Strings
    "'time' must be a single string" = is.character(time) && length(time) == 1,
    "'id' must be a single string or null" = is.null(id) || (is.character(id) && length(id) == 1),
    "'name' must be a character or null" = is.null(name) || is.character(name),

    # Consistency in given arguments
    "'nharm' cannot be greater than 'nbasis'" = nharm <= nbasis,
    "'dat' must have column names (the sites)" = !is.null(colnames(dat))
  )

  # ================================ TIEMPOS ===================================
  # Con la columna 'time' cada fila trae su tiempo, así que se pueden ordenar
  # sin riesgo. Sin ella el único orden posible es el de las filas.
  dat <- as.data.frame(dat, check.names = F)

  if (time %in% colnames(dat)) {
    times <- dat[[time]]
    dat   <- dat[, colnames(dat) != time, drop = F]
  } else {
    times <- seq_len(nrow(dat))
    message("No column '", time, "' in 'dat': the rows are taken as equally spaced times 1, ..., ",
            nrow(dat), ", in their current order.")
  }

  stopifnot(
    "The times must be numeric, Date or POSIXct" = is.numeric(times) || inherits(times, c('Date', 'POSIXt')),
    "The times cannot have missing values" = !anyNA(times),
    "The times must be unique" = !any(duplicated(times)),
    "'dat' must have at least one site besides the column of times" = ncol(dat) > 0
  )

  non_numeric <- colnames(dat)[!vapply(dat, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("The columns of the sites must be numeric. Not numeric: ", paste(non_numeric, collapse = ', '), call. = F)
  }

  ord   <- order(times)
  times <- times[ord]
  dat   <- as.matrix(dat[ord, , drop = F])
  rownames(dat) <- NULL
  t     <- as_days(times)

  # ============================== COORDENADAS =================================
  coords <- as.data.frame(coords, check.names = F)
  if (!is.null(id)) {
    if (!(id %in% colnames(coords))) stop("No column '", id, "' in 'coords'", call. = F)
    rownames(coords) <- as.character(coords[[id]])
    coords <- coords[, colnames(coords) != id, drop = F]
  }

  if (!all(vapply(coords, is.numeric, logical(1)))) {
    stop("The coordinates must be numeric", call. = F)
  }

  unknown <- setdiff(colnames(dat), rownames(coords))
  if (length(unknown) > 0) {
    stop("Sites of 'dat' without coordinates: ", paste(unknown, collapse = ', '), call. = F)
  }

  # SITIOS REPETIDOS: Dos estaciones en el mismo lugar dan dos filas iguales en
  # el sistema de kriging, que entonces no tiene solución
  site_coords <- as.matrix(coords[colnames(dat), , drop = F])
  repeated    <- duplicated(site_coords) | duplicated(site_coords, fromLast = T)
  if (any(repeated)) {
    groups <- split(colnames(dat)[repeated], apply(site_coords[repeated, , drop = F], 1, paste, collapse = ' '))
    stop("Sites with the same coordinates: ",
         paste(vapply(groups, paste, character(1), collapse = ' and '), collapse = '; '), call. = F)
  }

  # =========================== DATOS FALTANTES ================================
  # Se usan todos los datos disponibles: solo se quitan las estaciones sin
  # ninguna observación y las filas sin ninguna observación en ninguna estación.
  empty_sites <- colnames(dat)[colSums(!is.na(dat)) == 0]
  if (length(empty_sites) > 0) {
    message("Sites without observations, removed: ", paste(empty_sites, collapse = ', '))
    dat <- dat[, !(colnames(dat) %in% empty_sites), drop = F]
  }
  if (ncol(dat) == 0) stop("No site has observations", call. = F)

  observed_rows <- rowSums(!is.na(dat)) > 0
  dat   <- dat[observed_rows, , drop = F]
  times <- times[observed_rows]
  t     <- t[observed_rows]

  missing <- 100 * colMeans(is.na(dat))
  if (any(missing > 10)) {
    over <- missing[missing > 10]
    warning("Sites with more than 10% of missing values: ",
            paste0(names(over), ' (', round(over, 1), '%)', collapse = ', '), call. = F)
  }

  # Sin penalización una curva con menos observaciones que bases no tiene
  # solución única
  if (lambda == 0 && any(colSums(!is.na(dat)) < nbasis)) {
    stop("With 'lambda' = 0 every site needs at least 'nbasis' observations. ",
         "Use 'lambda' > 0 or fewer basis functions", call. = F)
  }

  # ========================== SUAVIZAMIENTO ===================================
  # Penalización de la segunda derivada (Ramsay y Silverman, 2005)
  hr     <- range(t)
  fbasis <- switch(basis,
                   'bsplines' = fda::create.bspline.basis(hr, nbasis),
                   'fourier'  = fda::create.fourier.basis(hr, nbasis),
                   stop("'basis' must be either 'bsplines' or 'fourier'", call. = F))
  fdpar  <- fda::fdPar(fdobj = fbasis, Lfdobj = fda::int2Lfd(2), lambda = lambda)

  fd <- smooth_curves(t, dat, fdpar)

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
  coords <- as.matrix(coords[colnames(dat), , drop = F])

  scores <- pca$scores
  dimnames(scores) <- list(colnames(dat), paste0(name, '.PC', seq_len(ncol(scores))))
  scores <- cbind(coords, scores)

  variance <- data.frame('basis'    = pca$harmonics$basis$type,
                         'no. PC'   = length(pca$varprop),
                         'variance' = sum(pca$varprop), check.names = F)
  rownames(variance) <- name

  # Lo funcional queda oculto, ver utils.R
  internal <- list(coords = coords, data = dat, time = times, t = t, fd = fd, pca = pca,
                   args = list(basis = basis, nbasis = nbasis, lambda = lambda))

  new_spatfd(scores, variance, missing, internal, name)
}
