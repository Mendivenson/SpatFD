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
#' The functional representation of the data is not needed for the spatial
#' modelling of the scores, but it is required later for prediction. It is
#' stored in environments attached as attributes, so that it is shared by
#' reference between `SpatFD` objects instead of being copied.
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
#' and the attributes (each one an environment)
#' \describe{
#'   \item{`pca`}{Output of [fda::pca.fd()].}
#'   \item{`fd`}{Functional data object, the `fd` element of
#'     [fda::smooth.basis()].}
#'   \item{`data`}{Original data matrix, plus the coordinates under `coords`.}
#'   \item{`models`}{Compiled variogram models, indexed by model name. See
#'     [add_model()].}
#' }
#' The first three are indexed by variable name.
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
spatfd <- function(dat,                      # Matriz/Dataframe de datos (Nombres columnas = estaciones)
                   coords,                   # Matriz coordendas (Nombres filas = estaciones)
                   t = 1:nrow(dat),          # Tiempos en que se hicieron las mediciones (3ra dimensión). Igual al número de filas del dataframe de dat
                   basis = 'bsplines',       # Base para el suavizamiento
                   nbasis = 91,              # Número de bases para el suavizamiento
                   lambda = 0.0002,          # Parámetro de suavizamiento
                   nharm = 20,               # Cantidad máxima de PC
                   var = NULL,               # Varianza deseada
                   name = NULL){             # Nombre de la variable

  # =================== SIMPLIFIED DATA VERIFICATION ===========================
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

  # ========================== DATA SMOOTHING ==================================
  hr <- range(t)                                          # Function's range
  oplfd <- fda::vec2Lfd(bwtvec = c(1,ncol(dat)),
                        rangeval = hr)                    # Penalization
  basis <- switch(basis,                                  # Functional basis to use
                  'bsplines' = fda::create.bspline.basis(hr,nbasis),
                  'fourier' = fda::create.fourier.basis(hr,nbasis),
                  stop("'basis' must be either 'bsplines' or 'fourier'", call. = F)
  )

  fun_par  <- fda::fdPar(fdobj=basis,Lfdobj=oplfd,lambda)
  smooth <- fda::smooth.basis(argvals=t,
                              dat,
                              fun_par)
  fd <- smooth$fd                                         # Functional data

  # ================================ PCA =======================================

  pca <- fda::pca.fd(fd,nharm=nharm)

  # AUTO PC SELECTION: If 'var' is specified, automatically determines the minimum
  # number of principal components required to meet the desired variance threshold.

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

  # --------------------------- RESULTS OBJECT ---------------------------------
  # When 'name' is NULL, the functional variable is named after the 'dat' argument using
  # deparsed substitution from the parent environment.

  if (is.null(name)){
    name <- deparse(substitute(dat))
  }

  result <- list()

  # Coordinates + scores resultantes del pca
  result$scores           <- pca$scores
  rownames(result$scores) <- colnames(dat)
  colnames(result$scores) <- paste0(name, '.PC', 1:ncol(pca$scores))
  result$scores           <- cbind(coords[colnames(dat),], result$scores)

  # Explained variance by the principal components
  result$variance <- data.frame('basis' = pca$harmonics$basis$type,
                                'no. PC' = length(pca$varprop),
                                'variance' = sum(pca$varprop), check.names = F)
  rownames(result$variance) <- name

  # Resumen de los modelos registrados: se calcula al consultarlo, ver `$.SpatFD`
  result['models'] <- list(NULL)

  # MEMORY OPTIMIZATION:
  # Although spatial analysis requires only the scores of the principal components, the underlying
  # functional representation (fd object) is necessary for prediction tasks. Objects not immediately
  # needed but required in subsequent stages are stored in separate environments. This design reduces
  # memory footprint by sharing references across multiple SpatFD objects rather than creating deep copies.

  env_pca <- new.env(parent = emptyenv())
  env_fd  <- new.env(parent = emptyenv())
  env_dat <- new.env(parent = emptyenv())

  env_pca[[name]]     <- pca
  env_fd[[name]]      <- fd
  env_dat[[name]]     <- dat
  env_dat[['coords']] <- coords

  result <- structure(result,
                      class = 'SpatFD',
                      pca = env_pca,
                      fd = env_fd,
                      data = env_dat,
                      models = new.env(parent = emptyenv()))

  return(result)
}

# ................................ SPATFD CLASS' METHODS ..........................................

# FUNCIONES DE IMPRESIÓN: Debido a la cartacterística del objeto, el usuario no
#                         puede acceder directamente a la forma funcional de los
#                         datos.
#
#                         summary.SpatFD() -> Imprime de una forma más organizada
#                                             dat$variance
#                         print.SpatFD()   -> Imprime dat$coords y dat$var_name*
#                                             de una forma "organizada"

#' Print and summarise a SpatFD object
#'
#' `print()` shows the coordinates and principal component scores of every
#' site. `summary()` shows, for each variable, the basis type, the number of
#' principal components and the proportion of variance they explain, followed
#' by the summary table of the registered models, if any.
#'
#' @param x,object An object of class `SpatFD`, as returned by [spatfd()].
#' @param ... Ignored.
#'
#' @return `print()` returns `x` invisibly. `summary()` is called for its side
#'   effect and returns `NULL` invisibly.
#'
#' @name SpatFD-print
#' @export
summary.SpatFD <- function(object, ...){
  cat('\n\n--------- PCA information ----------')
  knitr::kable(object$variance, format = 'pipe') |> print()
  cat('------------------------------------\n\n')

  if (nrow(object$models) > 0){
    cat('------------- Models ---------------')
    knitr::kable(object$models, format = 'pipe') |> print()
    cat('------------------------------------\n\n')
  }
}

# ACCESO A LOS CAMPOS: sfd$models no es un campo guardado sino el resumen, calculado
#                      al momento, de los modelos registrados en attr(sfd, 'models').
#                      Así la tabla nunca queda desactualizada cuando add_model()
#                      agrega un modelo en el sitio.

#' @export
`$.SpatFD` <- function(x, name){
  if (name == 'models') model_table(x) else .subset2(x, name)
}

#' @export
`[[.SpatFD` <- function(x, i, ...){
  if (identical(i, 'models')) model_table(x) else .subset2(x, i, ...)
}

#' @rdname SpatFD-print
#' @export
print.SpatFD <- function(x, ...){
  cat('\n\n-------------------------------------------------\n')
  print(x$scores)
  cat('-------------------------------------------------\n\n')
  invisible(x)
}

# FUNCIÓN DE PLOTEO: La función plot.SpatFD() tiene dos modos (Funcional y PCA)
#                    que respectivamente accedden a el método plot ya implementado
#                    en el paquete fda para los objetos provenientes de fda::smooth.basis()
#                    y fda::pca.fd()

#' Plot a SpatFD object
#'
#' Plots one of the functional variables stored in a `SpatFD` object, either
#' as smoothed curves or through its functional principal components. The
#' plotting itself is delegated to the methods of the \pkg{fda} package.
#'
#' @param x An object of class `SpatFD`, as returned by [spatfd()].
#' @param variable Name of the variable to plot. Defaults to the first
#'   variable of the object, see [var_names()].
#' @param mode Either `"functional"` to plot the smoothed curves or `"pca"` to
#'   plot the principal components.
#' @param multiple Passed as `ask` to the \pkg{fda} plot method when
#'   `mode = "functional"`.
#' @param ... Further arguments passed to the \pkg{fda} plot method.
#'
#' @return Called for its side effect of drawing a plot.
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
#' plot(sfd)
#'
#' @export
plot.SpatFD <- function(x, variable = NULL, mode = 'functional', multiple = F, ...){

  l <- ls(attr(x, which = 'fd'))
  if (is.null(variable)) {
    variable = l[1]
    message('Plotting the variable ', variable)
  }

  if (!(variable %in% l)){
    stop('Variable ', variable, ' is not registered in the object ', deparse(substitute(x)), call. = F)
  }

  switch(mode,
         'functional' =  {
           plot(attr(x, 'fd')[[variable]], ask = multiple, ...)
         },
         'pca' = plot(attr(x, 'pca')[[variable]], ...),
         stop("'mode' must be the string 'functional' or the string 'pca'", call. = F))
}

# var_name* representa todos los nomrbes de variables dados en name en algún momento
# del flujo de trabajo :D

# NOMBRES DE VARIABLES GUARDADAS EN EL OBJETO:

#' Names of the functional variables stored in an object
#'
#' @param x An object of class `SpatFD`, as returned by [spatfd()].
#'
#' @return A character vector with the names of the variables, in alphabetical
#'   order.
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
#' var_names(sfd)
#'
#' @export
var_names <- function(x) {
  UseMethod("var_names")
}

#' @rdname var_names
#' @export
var_names.SpatFD <- function(x) {
  ls(attr(x, 'fd'))
}


# OPERATOR OVERLOAD: A partir de dos objetos se puede crear un nuevo objeto multivariado

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
`+.SpatFD` <- function(x1,x2){

  if (class(x2) != 'SpatFD'){
    stop("Addition is only possible between two SpatFD objects.")
  }

  if (any(sort(colnames(attr(x1, 'data')$coords)) != sort(colnames(attr(x2, 'data')$coords)))){
    stop("Non compatible coordinates (Likely different names or dimensions)")
  }

  names_x1 <- var_names(x1)
  names_x2 <- var_names(x2)

  common_names <- intersect(names_x1, names_x2)

  if (length(common_names) > 0){
    warning('There are repeated variables between the two SpatFD objects. ',
            'Ignoring the repeated variables of the object ', deparse(substitute(x2)))
    names_x2 <- setdiff(names_x2, common_names)
  }

  result <- list()

  # Mezclar los dos litados de coordenadas reportados
  coords.x1  <- attr(x1, 'data')$coords
  coords.x2  <- attr(x2, 'data')$coords
  new_places <- setdiff(rownames(coords.x2), rownames(coords.x1))
  coords     <- rbind(coords.x1, coords.x2[new_places, , drop = F])

  # Extraer los scores
  score.x1 <- x1$scores |>
    (\(x) x[,-which(colnames(x1$scores) %in% colnames(coords)), drop = F])() |>
    as.data.frame()
  score.x2 <- x2$scores |>
    (\(x) x[,-which(colnames(x2$scores) %in% colnames(coords)), drop = F])() |>
    as.data.frame()

  # Reportar en una sola tabla la varianza retenida
  result$variance <- rbind(x1$variance, x2$variance[names_x2, , drop = F])

  # Outer join/full join
  pts <- c(rownames(score.x1), rownames(score.x2)) |>
    unique() |>
    sort()
  scores <- cbind(coords[pts,], score.x1[pts, , drop = F], score.x2[pts, , drop = F]) |>
    as.matrix()
  result$scores <- scores

  pca <- new.env(parent = emptyenv())
  fd <- new.env(parent = emptyenv())
  dat <- new.env(parent = emptyenv())

  list2env(as.list(attr(x1, 'pca')), envir = pca)
  list2env(as.list(attr(x2, 'pca'))[names_x2], envir = pca)

  list2env(as.list(attr(x1, 'fd')), envir = fd)
  list2env(as.list(attr(x2, 'fd'))[names_x2], envir = fd)

  list2env(as.list(attr(x1, 'data')), envir = dat)
  list2env(as.list(attr(x2, 'data'))[names_x2], envir = dat)

  dat$coords <- coords

  # Modelos registrados: primero los de x1 y luego los de x2, en su orden de registro
  models <- new.env(parent = emptyenv())
  for (x in list(x1, x2)){
    entries <- as.list(attr(x, 'models'))
    entries <- entries[order(vapply(entries, function(m) m$index, numeric(1)))]

    for (model in names(entries)){
      new_name <- model
      if (exists(model, envir = models, inherits = FALSE)){
        new_name <- next_model_name(models)
        message("Model '", model, "' of the second object renamed to '", new_name, "'")
      }
      entries[[model]]$index <- length(ls(models)) + 1
      models[[new_name]]     <- entries[[model]]
    }
  }

  result['models'] <- list(NULL)

  result <- structure(result,
                      class = 'SpatFD',
                      pca = pca,
                      fd = fd,
                      data = dat,
                      models = models)


  return(result)
  }
