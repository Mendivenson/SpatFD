# ================================ spatfd() ====================================
# Wraps the core FDA workflow (smoothing + PCA) from the `fda` package into a
# single structured object of class 'SpatFD'. The resulting object contains:
#
#   FIELDS (accessible via $):
#     $coords          -> Station coordinates with unique identifiers.
#     $variance        -> PCA summary: basis type, number of PCs, explained variance.
#     $<variable_name> -> PCA scores matrix (stations x PCs).
#
#   ATTRIBUTES (accessible via attr()):
#     attr(., 'pca')  -> Output of fda::pca.fd(). Stores harmonics and variance proportions.
#     attr(., 'fd')   -> Functional data object. Output of fda::smooth.basis()$fd.
#     attr(., 'data') -> Original raw data matrix.
#
# NOTE: Attributes are stored in separate environments to avoid unintended copies
#       under R's copy-on-modify semantics when the object is passed between functions.

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
      pca <- fda::pca.fd(fd, nharm = selected_pc[1])
      message('\t\033[32mVariance of ', var * 100, '% reached with ',
              selected_pc[1], ' principal components.\033[0m', appendLF = T)
    } else {
      pca <- fda::pca.fd(fd, nharm = nharm)
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
                      data = env_dat)

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
summary.SpatFD <- function(x){
  cat('\n\n--------- PCA information ----------')
  knitr::kable(x$variance, format = 'pipe') |> print()
  cat('------------------------------------\n\n')
}

print.SpatFD <- function(x){
  cat('\n\n-------------------------------------------------\n')
  print(x$scores)
  cat('-------------------------------------------------\n\n')
  invisible(x)
}

# FUNCIÓN DE PLOTEO: La función plot.SpatFD() tiene dos modos (Funcional y PCA)
#                    que respectivamente accedden a el método plot ya implementado
#                    en el paquete fda para los objetos provenientes de fda::smooth.basis()
#                    y fda::pca.fd()

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

var_names <- function(x) {
  UseMethod("var_names")
}
var_names.SpatFD <- function(x) {
  ls(attr(x, 'fd'))
}


# OPERATOR OVERLOAD: A partir de dos objetos se puede crear un nuevo objeto multivariado
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
    (\(x) x[,-which(colnames(x1$scores) %in% colnames(coords))])() |>
    as.data.frame()
  score.x2 <- x2$scores |>
    (\(x) x[,-which(colnames(x2$scores) %in% colnames(coords))])() |>
    as.data.frame()

  # Reportar en una sola tabla la varianza retenida
  result$variance <- rbind(x1$variance, x2$variance[names_x2, , drop = F])

  # Outer join/full join
  pts <- c(rownames(score.x1), rownames(score.x2)) |>
    unique() |>
    sort()
  scores <- cbind(coords[pts,], score.x1[pts,], score.x2[pts,]) |>
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

  result <- structure(result,
                      class = 'SpatFD',
                      pca = pca,
                      fd = fd,
                      data = dat)


  return(result)
  }
