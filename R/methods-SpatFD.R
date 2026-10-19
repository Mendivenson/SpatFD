# ============================ MÉTODOS DE SpatFD ===============================
# Lo que el usuario ve del objeto: los scores, la varianza explicada y la tabla
# de modelos. La forma funcional se usa para graficar, pero no se imprime.

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
#' @return `print()` returns `x` invisibly. `summary()` returns an object of
#'   class `summary.SpatFD`: a list with the tables `variance` and `models`.
#'
#' @name SpatFD-print
#' @export
print.SpatFD <- function(x, ...){
  cat('\n\n-------------------------------------------------\n')
  print(x$scores)
  cat('-------------------------------------------------\n\n')
  invisible(x)
}

#' @rdname SpatFD-print
#' @export
summary.SpatFD <- function(object, ...){
  structure(list(variance = object$variance, models = object$models),
            class = 'summary.SpatFD')
}

#' @rdname SpatFD-print
#' @export
print.summary.SpatFD <- function(x, ...){
  cat('\n--------- PCA information ----------\n')
  print(x$variance)
  cat('------------------------------------\n')

  if (nrow(x$models) > 0){
    cat('\n------------- Models ---------------\n')
    print(x$models)
    cat('------------------------------------\n')
  }
  invisible(x)
}

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
#'   plot the principal components, all of them on one page.
#' @param multiple Only for `mode = "functional"`: if `TRUE`, one curve per
#'   page, asking for `<Enter>` between them.
#' @param ... Further arguments passed to the \pkg{fda} plot method.
#'
#' @details
#' No plot of the package asks for `<Enter>` before drawing, even if the
#' graphics device was left in that mode (`par(ask = TRUE)`); the setting of
#' the device is restored afterwards. The only exception is
#' `multiple = TRUE`.
#'
#' @return `x`, invisibly. Called for its side effect of drawing a plot.
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

  if (is.null(variable)) {
    variable <- var_names(x)[1]
    message('Plotting the variable ', variable)
  }

  if (!(variable %in% var_names(x))){
    stop('Variable ', variable, ' is not registered in the object ', deparse(substitute(x)), call. = F)
  }

  # El gráfico no pide <Enter>, aunque el dispositivo esté en modo 'ask'. Solo
  # lo pide si el usuario lo solicita con multiple = TRUE (una curva por página).
  old <- par(ask = FALSE)
  on.exit(par(old), add = TRUE)

  switch(mode,
         'functional' = plot(get_fd(x, variable), ask = multiple, ...),
         'pca'        = plot_components(get_pca(x, variable), ...),
         stop("'mode' must be the string 'functional' or the string 'pca'", call. = F))

  invisible(x)
}

# COMPONENTES PRINCIPALES EN UNA SOLA PÁGINA: fda::plot.pca.fd() dibuja una
# componente por página y pide <Enter> entre ellas. Aquí se dibujan todas
# juntas, en una grilla de hasta tres columnas, salvo que el usuario ya haya
# dividido la página con par(mfrow).
plot_components <- function(pca, ...){
  K <- ncol(pca$harmonics$coefs)

  if (K > 1 && all(par('mfrow') == 1)) {
    columns <- min(K, 3)
    old <- par(mfrow = c(ceiling(K / columns), columns))
    on.exit(par(old), add = TRUE)
  }

  for (k in seq_len(K)) plot(pca, harm = k, ...)
}

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
  stopifnot('function only available for SpatFD objects' = inherits(x, 'SpatFD'))
  sort(names(attr(x, 'internal')$variables))
}
