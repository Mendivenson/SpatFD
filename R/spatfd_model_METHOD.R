# ========================== MÉTODOS DE SpatFDmodel ============================
# Lo que el usuario ve del modelo: los variogramas, la varianza que explican
# las componentes usadas y el gráfico de los variogramas.

#' Print and summarise a SpatFDmodel object
#'
#' `print()` shows the variograms of the model. `summary()` also shows, for
#' each variable, the number of components in the model and the proportion of
#' variance they explain.
#'
#' @param x,object An object of class `SpatFDmodel`, as returned by
#'   [spatfd_model()].
#' @param ... Ignored.
#'
#' @return `print()` returns `x` invisibly. `summary()` returns an object of
#'   class `summary.SpatFDmodel`: a list with the tables `variance` and
#'   `variograms`.
#'
#' @name SpatFDmodel-print
#' @export
print.SpatFDmodel <- function(x, ...){
  cat('\n------------- Variograms -------------\n')
  print(x$variograms)
  cat('--------------------------------------\n\n')
  invisible(x)
}

#' @rdname SpatFDmodel-print
#' @export
summary.SpatFDmodel <- function(object, ...){
  structure(list(variance = object$variance, variograms = object$variograms),
            class = 'summary.SpatFDmodel')
}

#' @rdname SpatFDmodel-print
#' @export
print.summary.SpatFDmodel <- function(x, ...){
  cat('\n------------- Variables --------------\n')
  print(x$variance)
  cat('--------------------------------------\n')
  cat('\n------------- Variograms -------------\n')
  print(x$variograms)
  cat('--------------------------------------\n')
  invisible(x)
}

#' Plot the variograms of one or several SpatFDmodel objects
#'
#' Draws the empirical variograms of the scores as points and the variograms
#' of the models as lines, in a matrix of panels like the one of \pkg{gstat}:
#' direct variograms on the diagonal and cross-variograms below it. Several
#' models can be drawn together to compare them, one line per model.
#'
#' @param x An object of class `SpatFDmodel`, as returned by [spatfd_model()].
#' @param y,... More `SpatFDmodel` objects to draw with `x`, and further
#'   named arguments passed to [gstat::variogram()], such as `cutoff` or
#'   `width`. The models are labelled with the names of their arguments
#'   (`plot(a, other = b)`) or, without a name, with the expression given.
#'
#' @details
#' The panels cover every score of every model. The empirical variograms are
#' computed with [gstat::variogram()] on the scores kept in the models, so
#' they are the same points seen when the models were fitted from
#' [to_gstat()]. Every panel has a legend; a model without a variogram for
#' that pair of scores is listed as not modelled.
#'
#' Without \pkg{gstat} only the lines are drawn, up to a third of the diagonal
#' of the area of the sites (the default cutoff of \pkg{gstat}).
#'
#' @return `x`, invisibly. Called for its side effect of drawing a plot.
#'
#' @examples
#' set.seed(1)
#' dat <- sapply(1:10, function(i) rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1))
#' colnames(dat) <- paste0("s", 1:10)
#' coords <- cbind(x = runif(10), y = runif(10))
#' rownames(coords) <- colnames(dat)
#' sim <- spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "sim")
#'
#' exponential <- data.frame(var = c("sim.PC1", "sim.PC2"), model = "exponential",
#'                           rango = 0.5, sigma = c(2, 0.5))
#' spherical   <- data.frame(var = "sim.PC1", model = "spherical", rango = 0.8, sigma = 2)
#'
#' plot(spatfd_model(sim, model = exponential))
#' plot(spatfd_model(sim, model = exponential), spherical = spatfd_model(sim, model = spherical))
#'
#' @export
plot.SpatFDmodel <- function(x, y, ...){

  # MODELOS Y ARGUMENTOS: Los objetos SpatFDmodel se dibujan; los demás
  # argumentos (con nombre) van a gstat::variogram()
  dots  <- list(...)
  args  <- c(list(x), if (!missing(y)) list(y), dots)
  exprs <- c(list(substitute(x)), if (!missing(y)) list(substitute(y)), as.list(substitute(list(...)))[-1])
  given <- c('', if (!missing(y)) '', if (is.null(names(dots))) rep('', length(dots)) else names(dots))

  is_model <- vapply(args, inherits, logical(1), what = 'SpatFDmodel')
  if (any(!is_model & given == '')) {
    stop('The unnamed arguments must be SpatFDmodel objects', call. = F)
  }

  models <- args[is_model]
  labels <- ifelse(given == '', vapply(exprs, function(e) paste(deparse(e), collapse = ''), character(1)), given)
  labels <- make.unique(labels[is_model])
  vario  <- args[!is_model]
  names(vario) <- given[!is_model]

  data   <- models_data(models)
  fields <- unlist(lapply(data, function(sfd) paste0(get_name(sfd), '.PC', seq_len(ncol(get_pca(sfd)$scores)))),
                   use.names = F)
  q      <- length(fields)

  # PUNTOS: Variogramas empíricos con gstat; sin gstat solo se dibujan las líneas
  if (requireNamespace('gstat', quietly = TRUE)) {
    points_df <- do.call(empirical_variograms, c(list(data), vario))
    cutoff    <- max(points_df$dist)
  } else {
    message("Without the package 'gstat' only the variograms of the models are drawn")
    points_df <- data.frame(dist = numeric(0), gamma = numeric(0), first = character(0), second = character(0))
    coords    <- do.call(rbind, lapply(data, get_coords))
    cutoff    <- sqrt(sum(apply(coords, 2, function(v) diff(range(v)))^2)) / 3
  }

  # Como gstat, las curvas empiezan justo después de cero: el efecto pepita se
  # ve como el intercepto
  h <- seq(cutoff / 200, cutoff, length.out = 200)

  # Un color y un tipo de línea por modelo, para distinguir líneas que coinciden
  colors <- rep(c('blue', 'red', 'darkgreen', 'darkorange', 'purple', 'brown'), length.out = length(models))
  types  <- seq_along(models)

  # MATRIZ DE PANELES: Directos en la diagonal y cruzados debajo, como en gstat
  old <- par(mfrow = c(q, q), mar = c(2.2, 2.2, 1.6, 0.4), mgp = c(1.2, 0.3, 0), tcl = -0.2,
             cex.axis = 0.75, ask = FALSE)
  on.exit(par(old), add = TRUE)

  for (i in seq_len(q)) {
    for (j in seq_len(q)) {
      if (j > i) {
        plot.new()
        next
      }

      a   <- fields[i]
      b   <- fields[j]
      pts <- points_df[(points_df$first == a & points_df$second == b) |
                       (points_df$first == b & points_df$second == a), ]
      fit <- lapply(models, function(model) {
        f <- pair_variogram(model, a, b)
        if (is.null(f)) NULL else f(h)
      })
      modelled <- !vapply(fit, is.null, logical(1))

      plot(NA, xlim = c(0, cutoff), ylim = range(0, pts$gamma, unlist(fit), na.rm = T),
           xlab = if (i == q) 'distance' else '', ylab = if (j == 1) 'semivariance' else '',
           main = if (i == j) a else paste(b, a, sep = '.'), cex.main = 0.85, font.main = 1)
      if (nrow(pts) > 0) points(pts$dist, pts$gamma, pch = 19, cex = 0.6, col = 'black')
      for (k in which(modelled)) lines(h, fit[[k]], lwd = 1.5, col = colors[k], lty = types[k])

      legend('bottomright', bty = 'n', cex = 0.65,
             legend = ifelse(modelled, labels, paste0(labels, ': not modelled')),
             col    = ifelse(modelled, colors, 'gray60'),
             lty    = ifelse(modelled, types, 0),
             lwd    = 1.5)
    }
  }

  invisible(x)
}
