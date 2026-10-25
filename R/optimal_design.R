# =============================== optimal_design() =============================

#' Optimal spatial design for functional data
#'
#' Finds where to place new stations so that the curves are predicted as
#' precisely as possible at a set of target points, given a model of the
#' variable and the stations it already has.
#'
#' @param model An object of class `SpatFDmodel` with a single variable, as
#'   returned by [spatfd_model()]. Its stations are the existing stations, and
#'   its components and variograms are the ones of the design.
#' @param k Number of new stations.
#' @param targets Where the prediction should be precise: a numeric matrix or
#'   data frame with the coordinates of the target points, or a map, an `sf`
#'   or `sfc` object with polygons, on which a regular grid of about
#'   `target_size` points is placed.
#' @param candidates Where the new stations can be placed: a numeric matrix or
#'   data frame with their coordinates, or `"random"` (the default) for
#'   `candidate_size` random points, inside the map when `targets` is a map
#'   and inside the rectangle that contains the target points and the
#'   stations otherwise. Use [set.seed()] before to repeat the random points.
#'   Candidates at the place of an existing station are dropped.
#' @param method Prediction method whose variance is minimised: `"scores"`
#'   (the default) or `"lambda"`, see [predict()][predict.SpatFDmodel].
#' @param target_size Only when `targets` is a map: approximate number of
#'   points of its grid.
#' @param candidate_size Only when `candidates` is `"random"`: number of
#'   random candidates.
#' @param progress If `TRUE`, progress bars with the candidates tried in each
#'   stage: preparing the candidates, adding the stations and each round of
#'   exchanges. By default, only in interactive sessions.
#'
#' @details
#' **What is minimised.** The prediction variance of the curve,
#' \eqn{\sigma^2(s_0)}, added up over the target points:
#' \deqn{V = \sum_{s_0} \sigma^2(s_0)}
#' It is the variance returned by [predict()][predict.SpatFDmodel], with
#' ordinary kriging: with the method `"scores"`, the sum of the variances of
#' the scores; with `"lambda"`, the variance of the curve. The variance does
#' not depend on the observed values, only on where the stations are, so a
#' design can be evaluated before measuring.
#'
#' **How it is searched.** Among the candidates, in two steps:
#' \enumerate{
#'   \item Greedy: the stations are added one at a time, each time the
#'     candidate that most reduces \eqn{V}.
#'   \item Exchanges: each new station is swapped with every free candidate,
#'     keeping the swap when \eqn{V} decreases, until no swap improves it.
#' }
#' The search is deterministic, the stations are always candidates (inside
#' the map, never repeated) and the variance reported is the one of the final
#' design.
#'
#' **Speed.** Everything that does not change between designs (the system of
#' the existing stations and the variograms between candidates and targets)
#' is computed once; evaluating a design then only inverts a small matrix of
#' the size of the design. The time grows mostly with `candidate_size` and `k`;
#' the target points are cheap. With the defaults, around ten seconds.
#'
#' Only models of a single variable are supported for now.
#'
#' @return An object of class `SpatFDdesign`: a list with
#' \describe{
#'   \item{`new_stations`}{Matrix with the coordinates of the new stations,
#'     in the order they were added.}
#'   \item{`variance`}{Data frame with the total variance \eqn{V} and its mean
#'     over the target points: without new stations, after adding each one,
#'     and after the exchanges.}
#' }
#' It can be plotted with [plot()][plot.SpatFDdesign].
#'
#' @references
#' Bohórquez, M., Giraldo, R. and Mateu, J. (2016). Optimal sampling for
#' spatial prediction of functional data. *Statistical Methods &
#' Applications*, 25(1), 39-54. \doi{10.1007/s10260-015-0340-9}
#'
#' @seealso [spatfd_model()], [predict()][predict.SpatFDmodel]
#'
#' @examples
#' set.seed(1)
#' dat <- sapply(1:10, function(i) rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1))
#' colnames(dat) <- paste0("s", 1:10)
#' coords <- cbind(x = runif(10, 0, 0.5), y = runif(10))
#' rownames(coords) <- colnames(dat)
#' sim <- spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "sim")
#'
#' model <- spatfd_model(sim, model = data.frame(var = c("sim.PC1", "sim.PC2"), model = "exponential",
#'                                               rango = 0.3, sigma = c(30, 0.1)))
#'
#' # Targets on a regular grid over the unit square
#' targets <- expand.grid(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1))
#' design  <- optimal_design(model, k = 3, targets = targets, candidate_size = 200)
#' design
#' plot(design)
#'
#' @export
optimal_design <- function(model,                      # Modelo de una variable, ver spatfd_model()
                           k,                          # Número de estaciones nuevas
                           targets,                    # Puntos objetivo (coordenadas o mapa)
                           candidates = 'random',      # Candidatos (coordenadas o 'random')
                           method = 'scores',          # Método cuya varianza se minimiza
                           target_size = 2000,         # Puntos de la grilla, si 'targets' es un mapa
                           candidate_size = 5000,      # Candidatos al azar, si candidates = 'random'
                           progress = interactive()){  # Barra de progreso

  stopifnot(
    "'model' must be of class SpatFDmodel" = inherits(model, 'SpatFDmodel'),
    "'progress' must be TRUE or FALSE" = isTRUE(progress) || isFALSE(progress),
    "'k' must be a positive whole number" = is.numeric(k) && length(k) == 1 && k >= 1 && k == round(k),
    "'method' must be a single string" = is.character(method) && length(method) == 1,
    "'candidate_size' must be a positive whole number" =
      is.numeric(candidate_size) && length(candidate_size) == 1 && candidate_size >= 1
  )
  method <- match.arg(method, c('scores', 'lambda'))

  data <- get_model_data(model)
  if (length(data) > 1) {
    stop('optimal_design() only supports models of a single variable for now', call. = F)
  }
  stations    <- get_coords(data[[1]])
  coord_names <- colnames(stations)

  # PUNTOS OBJETIVO: Coordenadas, o una grilla regular sobre un mapa
  map <- NULL
  if (inherits(targets, c('sf', 'sfc'))) {
    map     <- targets
    targets <- map_grid(map, target_size)
  }
  targets <- check_newcoords(targets, coord_names)

  # CANDIDATOS: Coordenadas o al azar, sin los que caen sobre una estación
  if (identical(candidates, 'random')) {
    candidates <- random_candidates(candidate_size, map, targets, stations)
  }
  candidates <- check_newcoords(candidates, coord_names)
  on_station <- apply(candidates, 1, function(p) any(colSums(abs(t(stations) - p)) == 0))
  candidates <- unique(candidates[!on_station, , drop = F])

  if (k > nrow(candidates)) {
    stop("'k' cannot be greater than the number of candidates (", nrow(candidates), ')', call. = F)
  }

  # BÚSQUEDA: Voraz y después intercambios, con la varianza rápida de los
  # diseños (todo lo que no cambia entre candidatos se calcula una vez), ver
  # optimal_design_UTILS.R
  systems <- get_compiled(model)[[method]]
  if (method == 'lambda') systems <- systems$systems
  N       <- nrow(candidates)
  cache   <- design_cache(systems, candidates, targets, show = progress)

  initial  <- cache$value(integer(0))
  greedy   <- greedy_design(cache, k, N, show = progress)
  exchange <- exchange_design(cache, greedy$chosen, greedy$variance[k], N, show = progress)

  new_stations <- candidates[exchange$chosen, , drop = F]
  rownames(new_stations) <- paste0('new', seq_len(k))

  total    <- c(initial, greedy$variance, exchange$variance)
  variance <- data.frame(step     = c('existing stations', paste0('+ new', seq_len(k)), 'after exchanges'),
                         stations = c(nrow(stations), nrow(stations) + seq_len(k), nrow(stations) + k),
                         total    = total,
                         mean     = total / nrow(targets))

  structure(list(new_stations = new_stations, variance = variance),
            class      = 'SpatFDdesign',
            stations   = stations,
            targets    = targets,
            candidates = candidates,
            map        = map,
            method     = method,
            variable   = names(data))
}
