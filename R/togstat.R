#' Export the scores of a SpatFD object to gstat
#'
#' Builds a `gstat` object with one entry per principal component score of the
#' selected variables, each with a constant mean (`score ~ 1`). The result can
#' be used with [gstat::variogram()], [gstat::fit.variogram()] or
#' [gstat::fit.lmc()].
#'
#' @param sfd An object of class `SpatFD`, as returned by [spatfd()].
#' @param vars Character vector with the names of the variables to export.
#'   Defaults to all the variables in `sfd`.
#'
#' @details
#' Sites with a missing score are dropped from the corresponding entry, so
#' variables observed at different sets of sites can be exported together.
#'
#' @return An object of class `gstat` whose ids are the names of the score
#'   columns (`<variable>.PC1`, `<variable>.PC2`, ...).
#'
#' @section Theory:
#' After the functional principal component analysis, the spatial dependence
#' between curves is determined by the spatial dependence of the scores
#' \eqn{f_k(s)}, which are scalar random fields. Each score is therefore
#' treated as a classical geostatistical variable: one variogram per score in
#' the univariate case, plus the cross-variograms between scores when several
#' functional variables are modelled together (Bohorquez, Giraldo and Mateu,
#' 2016, 2017). \pkg{gstat} (Pebesma, 2004) is one tool to estimate and fit
#' them.
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
#' Pebesma, E. J. (2004). Multivariable geostatistics in S: the gstat package.
#' *Computers & Geosciences*, 30, 683-691. \doi{10.1016/j.cageo.2004.03.012}
#'
#' @seealso [toSpatFD()] to translate the fitted models back.
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
#' g <- togstat(sfd)
#' names(g$data)
#'
#' @importFrom stats reformulate
#' @export
togstat <- function(sfd, vars = var_names(sfd)) {

  # CLASS HANDLING
  stopifnot(
    'function only available for SpatFD objects' = inherits(sfd, 'SpatFD'),
    'Variables not found in object' = all(vars %in% var_names(sfd))
  )

  # SEPARAR SCORES Y COORDENADAS:
  # Se separan de esta forma para permitir cualquier nombre de coordenada y variable posible
  scores <- sfd$scores[, !(colnames(sfd$scores) %in% colnames(attr(sfd, 'data')$coords))]
  coords <- attr(sfd, 'data')$coords

  # SELECCIÓN DE LOS COMPONENTES PRINCIPALES CON LOS QUE SE VA A TRABAJAR:
  # Básicamente se seleccionan las variables escogidas por el usuario.
  scores <- scores[, grepl(paste0("^(", paste0(vars, '\\.', collapse = '|'), ")"),
                           colnames(scores)), drop = FALSE]

  # SE AÑADEN LAS COORDENADAS ESPACIALES A LOS SCORES.
  scores <- as.data.frame(scores)
  sp::coordinates(scores) <- coords[rownames(scores), ]

  # GENERACIÓN OBJETO GSTAT:
  # El objeto gstat se genera en conjunto independientemente de si se seleccionan una o varias variables
  g <- NULL
  for (j in names(scores)) {
    g <- gstat::gstat(g       = g,
                      id      = j,
                      formula = reformulate('1', response = j),
                      data    = scores[!is.na(scores[[j]]), j])
  }

  return(g)
}
