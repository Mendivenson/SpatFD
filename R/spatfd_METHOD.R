# ============================ MÉTODOS DE SpatFD ===============================
# Lo que el usuario ve del objeto: los scores, la varianza explicada y los
# datos faltantes. La forma funcional se usa para graficar, pero no se imprime.

#' Print and summarise a SpatFD object
#'
#' `print()` shows the coordinates and principal component scores of every
#' site. `summary()` shows the basis type, the number of principal components
#' and the proportion of variance they explain, followed by the percentage of
#' missing values of the sites that have any.
#'
#' @param x,object An object of class `SpatFD`, as returned by [spatfd()].
#' @param ... Ignored.
#'
#' @return `print()` returns `x` invisibly. `summary()` returns an object of
#'   class `summary.SpatFD`: a list with the elements `variance` and
#'   `missing`.
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
  structure(list(variance = object$variance, missing = object$missing),
            class = 'summary.SpatFD')
}

#' @rdname SpatFD-print
#' @export
print.summary.SpatFD <- function(x, ...){
  cat('\n--------- PCA information ----------\n')
  print(x$variance)
  cat('------------------------------------\n')

  # Solo se listan los sitios con faltantes
  missing <- x$missing[x$missing > 0]
  if (length(missing) > 0){
    cat('\n------ Missing values (%) ----------\n')
    print(round(missing, 2))
    cat('------------------------------------\n')
  }
  invisible(x)
}

#' Plot a SpatFD object
#'
#' Plots the functional variable of a `SpatFD` object, either as smoothed
#' curves or through its functional principal components. The plotting itself
#' is delegated to the methods of the \pkg{fda} package.
#'
#' @param x An object of class `SpatFD`, as returned by [spatfd()].
#' @param mode Either `"functional"` to plot the smoothed curves or `"pca"` to
#'   plot the principal components, all of them on one page.
#' @param multiple Only for `mode = "functional"`: if `TRUE`, one curve per
#'   page, asking for `<Enter>` between them.
#' @param ... Further arguments passed to the \pkg{fda} plot method.
#'
#' @details
#' The horizontal axis is the time used in the smoothing: the times of `dat`
#' when they are numbers, and days since 1970-01-01 when they are dates, see
#' [spatfd()].
#'
#' No plot of the package asks for `<Enter>` before drawing, even if the
#' graphics device was left in that mode (`par(ask = TRUE)`); the setting of
#' the device is restored afterwards. The only exception is
#' `multiple = TRUE`.
#'
#' @return `x`, invisibly. Called for its side effect of drawing a plot.
#'
#' @examples
#' set.seed(1)
#' dat <- sapply(1:10, function(i) rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1))
#' colnames(dat) <- paste0("s", 1:10)
#' coords <- cbind(x = runif(10), y = runif(10))
#' rownames(coords) <- colnames(dat)
#' sfd <- spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "sim")
#'
#' plot(sfd)
#' plot(sfd, mode = "pca")
#'
#' @export
plot.SpatFD <- function(x, mode = 'functional', multiple = F, ...){

  # El gráfico no pide <Enter>, aunque el dispositivo esté en modo 'ask'. Solo
  # lo pide si el usuario lo solicita con multiple = TRUE (una curva por página).
  old <- par(ask = FALSE)
  on.exit(par(old), add = TRUE)

  switch(mode,
         'functional' = plot(get_fd(x), ask = multiple, ...),
         'pca'        = plot_components(get_pca(x), ...),
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
