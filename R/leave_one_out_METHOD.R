# =========================== MÉTODOS DE SpatFDloo =============================

#' Print a SpatFDloo object
#'
#' Shows, for each variable, the table of the cross-validation errors: one row
#' per model and method, from the smallest to the largest mean error.
#'
#' @param x An object of class `SpatFDloo`, as returned by [leave_one_out()].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.SpatFDloo <- function(x, ...){
  cat('\n------ Leave-one-out cross-validation ------\n')
  cat('Root mean square error over time, in the units of the data\n')
  for (variable in names(x$summary)) {
    cat('\n', variable, ':\n', sep = '')
    print(x$summary[[variable]], digits = 4, row.names = FALSE)
  }
  cat('--------------------------------------------\n\n')
  invisible(x)
}

#' Plot the error curves of a cross-validation
#'
#' Draws, for one variable, the error curves of the cross-validation in a
#' grid of panels, one per model and method: in each panel, how the error
#' \eqn{e(t) = \chi(t) - \hat\chi(t)} of every site evolves in time when the
#' site was left out.
#'
#' @param x An object of class `SpatFDloo`, as returned by [leave_one_out()].
#' @param variable Name of the variable to plot. Defaults to the first
#'   validated variable.
#' @param ... Graphical parameters passed to the lines, for example `lwd`.
#'
#' @details
#' Each site has the same colour in every panel, and all the panels share the
#' vertical scale, so the models and methods can be compared at a glance. A
#' dashed line marks the zero error. When the times of the data are dates,
#' the horizontal axis shows them as dates.
#'
#' @return `x`, invisibly. Called for its side effect of drawing a plot.
#'
#' @seealso [leave_one_out()]
#'
#' @examples
#' set.seed(1)
#' dat <- sapply(1:10, function(i) rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1))
#' colnames(dat) <- paste0("s", 1:10)
#' coords <- cbind(x = runif(10), y = runif(10))
#' rownames(coords) <- colnames(dat)
#' sim <- spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "sim")
#'
#' model <- spatfd_model(sim, model = data.frame(var = c("sim.PC1", "sim.PC2"), model = "exponential",
#'                                               rango = 0.5, sigma = c(30, 0.1)))
#' plot(leave_one_out(model))
#'
#' @export
plot.SpatFDloo <- function(x, variable = NULL, ...){

  internal <- attr(x, 'internal')
  if (is.null(variable)) {
    variable <- names(internal$curves)[1]
    if (length(internal$curves) > 1) message('Plotting the variable ', variable)
  }
  if (!(variable %in% names(internal$curves))) {
    stop('Variable ', variable, ' was not validated. The validated variables are: ',
         paste(names(internal$curves), collapse = ', '), call. = F)
  }

  curves   <- internal$curves[[variable]]
  times    <- internal$times[[variable]]
  template <- internal$template[[variable]]

  # COLORES: Uno por sitio, el mismo en todos los paneles
  sites  <- unique(unlist(lapply(curves, colnames)))
  colors <- stats::setNames(hcl.colors(length(sites), 'Dark 3'), sites)

  # GRILLA DE PANELES: Uno por modelo y método, con la misma escala vertical
  panels  <- names(curves)
  columns <- min(length(panels), 2)
  limits  <- range(unlist(curves))

  old <- par(mfrow = c(ceiling(length(panels) / columns), columns), mar = c(3.5, 3.5, 2, 0.5),
             mgp = c(2, 0.6, 0), ask = FALSE)
  on.exit(par(old), add = TRUE)

  style <- utils::modifyList(list(lty = 1, lwd = 1.2), list(...))

  for (panel in panels) {
    # El eje se arma con los tiempos del usuario; las líneas usan su valor
    # numérico interno (días para Date, segundos para POSIXct)
    axis_t <- from_days(times[[panel]], template)

    plot(range(axis_t), limits, type = 'n', xlab = 'time', ylab = 'error',
         main = sub('\\.(scores|lambda)$', ' - \\1', panel), font.main = 1, cex.main = 1)
    abline(h = 0, lty = 2, col = 'grey50')
    do.call(matlines, c(list(as.numeric(axis_t), curves[[panel]], col = colors[colnames(curves[[panel]])]), style))
  }

  invisible(x)
}
