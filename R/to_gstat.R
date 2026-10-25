# ================================= to_gstat() =================================

#' Export the scores of SpatFD objects to gstat
#'
#' Builds a `gstat` object with one entry per principal component score of the
#' given variables, each with a constant mean (`score ~ 1`), ready to estimate
#' and fit variograms with [gstat::variogram()], [gstat::fit.variogram()] or
#' [gstat::fit.lmc()]. It is a convenience for the user: the fitted model is
#' then given to [spatfd_model()] together with the same `SpatFD` objects.
#'
#' @param ... Objects of class `SpatFD`, as returned by [spatfd()], with
#'   different names (the attribute `name`).
#' @param pcs Optional named vector with the number of principal components to
#'   export for each variable, named after the variables: the first `pcs`
#'   components are exported. Variables not in `pcs` export all their
#'   components.
#'
#' @details
#' Each variable keeps its own sites, so variables observed at different
#' stations can be exported together. The entries follow the order of the
#' objects in `...` and, within each variable, the order of the components;
#' [gstat::fit.lmc()] depends on that order.
#'
#' The coordinates of all the objects must have the same column names.
#' Requires the packages \pkg{gstat} and \pkg{sp}.
#'
#' @return An object of class `gstat` whose ids are the names of the scores,
#'   `<name>.PC1`, `<name>.PC2`, ...
#'
#' @section Theory:
#' After the functional principal component analysis, the spatial dependence
#' between curves is carried by the scores \eqn{f_k(s)}, which are scalar
#' random fields. Each score is then a classical geostatistical variable: it
#' gets its own variogram, plus cross-variograms with other scores when
#' several of them are modelled jointly (Bohórquez, Giraldo and Mateu, 2016,
#' 2017).
#'
#' @references
#' Bohórquez, M., Giraldo, R. and Mateu, J. (2016). Optimal sampling for
#' spatial prediction of functional data. *Statistical Methods &
#' Applications*, 25(1), 39-54. \doi{10.1007/s10260-015-0340-9}
#'
#' Bohórquez, M., Giraldo, R. and Mateu, J. (2017). Multivariate functional
#' random fields: prediction and optimal sampling. *Stochastic Environmental
#' Research and Risk Assessment*, 31, 53-70. \doi{10.1007/s00477-016-1266-y}
#'
#' Pebesma, E. J. (2004). Multivariable geostatistics in S: the gstat package.
#' *Computers & Geosciences*, 30, 683-691. \doi{10.1016/j.cageo.2004.03.012}
#'
#' @seealso [spatfd_model()]
#'
#' @examplesIf requireNamespace("gstat", quietly = TRUE)
#' data(air_mexico)
#' pm10 <- spatfd(pm10_mexico, coords_pm10_mexico, basis = "fourier",
#'                nbasis = 21, lambda = 0.000001, nharm = 3, name = "PM10")
#' no2  <- spatfd(no2_mexico, coords_no2_mexico, basis = "fourier",
#'                nbasis = 21, lambda = 0.000001, nharm = 3, name = "NO2")
#'
#' g <- to_gstat(pm10, no2, pcs = c(PM10 = 2, NO2 = 1))
#' names(g$data)
#'
#' @export
to_gstat <- function(...,                    # Objetos SpatFD
                     pcs = NULL){            # Número de PC por variable (con nombres)

  if (!requireNamespace('gstat', quietly = TRUE) || !requireNamespace('sp', quietly = TRUE)) {
    stop("to_gstat() needs the packages 'gstat' and 'sp'", call. = F)
  }

  objects <- spatfd_list(...)
  pcs     <- check_pcs(objects, pcs)

  # Una entrada por score, cada una con los sitios de su variable
  g <- NULL
  for (name in names(objects)) {
    sfd    <- objects[[name]]
    coords <- get_coords(sfd)

    for (k in seq_len(pcs[[name]])) {
      id <- paste0(name, '.PC', k)
      scores <- data.frame(coords, get_pca(sfd)$scores[, k], check.names = F)
      colnames(scores)[ncol(scores)] <- id
      sp::coordinates(scores) <- colnames(coords)

      g <- gstat::gstat(g       = g,
                        id      = id,
                        formula = stats::as.formula(call('~', as.name(id), 1)),
                        data    = scores)
    }
  }

  g
}
