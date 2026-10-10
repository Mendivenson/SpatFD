#' Variogram models available in SpatFD
#'
#' Named list of variogram model functions. It mirrors the models of
#' \pkg{gstat} and adds some that \pkg{gstat} does not provide.
#'
#' @format A named list of functions. Every function takes a numeric vector of
#' distances `h` as its first argument and returns the variogram at those
#' distances. Depending on the model, it also takes
#' \describe{
#'   \item{`sigma`}{Partial sill (the slope or scale for unbounded models).}
#'   \item{`rango`}{Range parameter.}
#'   \item{`alpha1`}{Additional shape parameter (`kappa` in \pkg{gstat}).}
#' }
#'
#' Every model checks its own parameters and stops with an error when one is
#' outside its domain: `h` must be numeric, non-negative and without missing values;
#' `sigma` a single finite number; `rango` and `alpha1` single numbers greater
#' than zero. In addition, `rango` (the exponent) of `power` and `alpha1` of
#' `stable` cannot exceed 2, and `rango` of `linear` can be zero, which gives
#' the linear model without a sill. The sign of `sigma` is not checked, because
#' the cross-variograms of a linear model of coregionalisation can have
#' negative partial sills.
#'
#' Mirrors of \pkg{gstat} models (\pkg{gstat} code in parentheses):
#' `nugget` (`Nug`), `exponential` (`Exp`), `spherical` (`Sph`), `gaussian`
#' (`Gau`), `stable` (`Exc`), `matern` (`Mat`), `stein's matern` (`Ste`),
#' `circular` (`Cir`), `linear` (`Lin`), `bessel` (`Bes`), `pentamodel`
#' (`Pen`), `periodic` (`Per`), `wave` (`Wav`), `hole` (`Hol`), `logarithmic`
#' (`Log`), `power` (`Pow`), `spline` (`Spl`) and `legendre` (`Leg`). They
#' follow the formulas of the source code of \pkg{gstat} and reproduce its
#' values. Two of them have peculiarities inherited from it: `spline` returns
#' the sill beyond its range, and `legendre` takes the distances in kilometres
#' over the sphere.
#'
#' Models not available in \pkg{gstat}: `lanteujoul`, `radon ord 2`,
#' `radon ord 4` and `cubic`.
#'
#' @section Theory:
#' The variogram of a score field \eqn{f(s)} measures how different its values
#' are at two sites separated by a distance \eqn{h}:
#' \deqn{\gamma(h) = \frac{1}{2} Var\left( f(s + h) - f(s) \right)}
#' A variogram model is a function of \eqn{h} that is guaranteed to be valid
#' (conditionally negative definite). For bounded models, `sigma` is the value
#' at which the variogram levels off and `rango` controls the distance at
#' which it does so; beyond it, observations are practically uncorrelated
#' (Chilès and Delfiner, 2012). The mirrored models are meant to follow the
#' parametrisation of \pkg{gstat} (Pebesma, 2004).
#'
#' @references
#' Chilès, J.-P. and Delfiner, P. (2012). *Geostatistics: Modeling Spatial
#' Uncertainty* (2nd ed.). John Wiley & Sons. \doi{10.1002/9781118136188}
#'
#' Pebesma, E. J. (2004). Multivariable geostatistics in S: the gstat package.
#' *Computers & Geosciences*, 30, 683-691. \doi{10.1016/j.cageo.2004.03.012}
#'
#' @seealso [availablemodels()] to list or plot the models, and
#'   [add_model()] to register models built from them in a `SpatFD` object.
#'
#' @examples
#' names(models)
#' models$exponential(h = c(0, 1, 5), rango = 1.5, sigma = 5)
#'
#' @export
models <- list(

  # ----------------------------------- GSTAT SUPPORTED MODELS -------------------------------------

  # ::efecto pepita::
  nugget = function(h, sigma = 2, ...){
    stopifnot(
      'Model "nugget": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "nugget": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    sigma * (h > 0)        # El variograma vale cero en el origen
  },

  # ::modelo exponencial::
  exponential = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "exponential": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "exponential": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "exponential": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    sigma - sigma * exp(- h/rango)
  },

  # ::modelo esférico::
  spherical = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "spherical": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "spherical": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "spherical": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    vg     <- h
    bfr    <- vg < rango                   # Antes de alcanzar el rango
    afr    <- !bfr                         # Después de alcanzar el rango

    vg[bfr] <- sigma - sigma * (1 - 1.5 * h[bfr]/rango + 0.5 * (h[bfr]/rango)^3)
    vg[afr] <- sigma
    vg
  },

  # ::modelo gaussiano::
  gaussian = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "gaussian": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "gaussian": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "gaussian": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    sigma - sigma * exp(- ((h/rango)^2))
  },

  # ::modelo estable::
  stable = function(h, rango = 1, sigma = 2, alpha1 = 1.3,...){
    stopifnot(
      'Model "stable": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "stable": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "stable": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma),
      'Model "stable": `alpha1` must be a single number in (0, 2]' =
        is.numeric(alpha1) && length(alpha1) == 1 && is.finite(alpha1) && alpha1 > 0 && alpha1 <= 2
    )
    vg  <- h
    vg  <- sigma - sigma * exp( -((h/rango)^alpha1))
    vg
  },

  # ::modelo matern::
  matern = function(h, rango = 1, sigma = 2, alpha1 = 0.25, ...){
    stopifnot(
      'Model "matern": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "matern": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "matern": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma),
      'Model "matern": `alpha1` must be a single number greater than 0' =
        is.numeric(alpha1) && length(alpha1) == 1 && is.finite(alpha1) && alpha1 > 0
    )
    vg <- sigma - sigma *
      (1 / ((2**(alpha1-1)) * gamma(alpha1))) * ((h/rango)**alpha1) *
      besselK(x = h/rango, nu = alpha1)
    vg[h == 0] <- 0         # EVITAR: Nan en Bessel eval
    vg
  },

  # ::modelo matern (parametrización Stein)::
  `stein's matern` = function(h, rango = 1, sigma = 2, alpha1 = 0.25, ...){
    stopifnot(
      "Model \"stein's matern\": `h` must be numeric, non-negative and without missing values" =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      "Model \"stein's matern\": `rango` must be a single number greater than 0" =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      "Model \"stein's matern\": `sigma` must be a single finite number" =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma),
      "Model \"stein's matern\": `alpha1` must be a single number greater than 0" =
        is.numeric(alpha1) && length(alpha1) == 1 && is.finite(alpha1) && alpha1 > 0
    )
    arg <- 2 * sqrt(alpha1) * h / rango
    vg <- sigma - sigma *
      (2^(1 - alpha1) / gamma(alpha1)) * (arg^alpha1) *
      besselK(x = arg, nu = alpha1)
    vg[h == 0] <- 0        # EVITAR: Nan en Bessel eval
    vg
  },

  # ::modelo circular::
  circular = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "circular": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "circular": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "circular": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    vg  <- h
    bfr <- vg < rango       # Antes de alcanzar el rango
    afr <- !bfr             # Después de alcanzar el rango

    hr      <- h[bfr]/rango
    vg[bfr] <- sigma * (2/pi) * (hr * sqrt(1 - hr^2) + asin(hr))
    vg[afr] <- sigma
    vg
  },

  # ::modelo lineal::
  # Con rango = 0 es el modelo lineal sin meseta; con rango > 0 crece
  # linealmente hasta alcanzar la meseta en el rango
  linear = function(h, rango = 0, sigma = 2, ...){
    stopifnot(
      'Model "linear": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "linear": `rango` must be a single number greater than or equal to 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango >= 0,
      'Model "linear": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    if (rango == 0) sigma * h else sigma * pmin(h/rango, 1)
  },

  # ::modelo Bessel::
  bessel = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "bessel": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "bessel": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "bessel": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    vg <- sigma - sigma * (h/rango) * besselK(x = h/rango, nu = 1)
    vg[h == 0] <- 0         # EVITAR: Nan en Bessel eval
    vg
  },

  # ::modelo pentaesférico::
  pentamodel = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "pentamodel": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "pentamodel": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "pentamodel": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    vg  <- h
    bfr <- vg < rango          # Antes de alcanzar el rango
    afr <- !bfr                # Luego de alcanzar el rango

    hr      <- h[bfr]/rango
    vg[bfr] <- sigma * hr * (15/8 - (5/4) * hr^2 + (3/8) * hr^4)
    vg[afr] <- sigma
    vg
  },

  # ::modelo periódico::
  periodic = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "periodic": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "periodic": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "periodic": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    sigma - sigma * cos((2 * pi * h) / rango)
  },

  # ::modelo wave::
  wave        = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "wave": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "wave": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "wave": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    vg  <- h

    vg  <- sigma - sigma * (rango/ (pi * h)) * sin(pi * h/rango)
    vg[h == 0] <- 0
    vg
  },

  # ::modelo hueco::
  hole = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "hole": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "hole": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "hole": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    vg <- sigma - sigma * sin(h/rango) / (h/rango)
    vg[h == 0] <- 0    # EVITAR: 0/0
    vg
  },

  # ::modelo logarítmico::
  logarithmic = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "logarithmic": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "logarithmic": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "logarithmic": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    vg <- sigma * log(h + rango)
    vg[h == 0] <- 0    # El variograma vale cero en el origen
    vg
  },

  # ::modelo power::
  power = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "power": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "power": `rango` must be a single number in (0, 2]' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0 && rango <= 2,
      'Model "power": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    sigma * (h ^ rango)
  },

  # ::modelo spline::
  # gstat trata este modelo como transitivo: más allá del rango devuelve la
  # meseta. En un modelo anidado gstat aplica ese corte en el mayor de los
  # rangos de todas las estructuras (y no lo aplica si alguna alcanza la meseta
  # de forma asintótica); aquí el corte es siempre en el rango de esta estructura.
  spline = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "spline": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "spline": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "spline": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    vg <- sigma * h * h * log(h)
    vg[h == 0]    <- 0        # EVITAR: -Inf log eval
    vg[h > rango] <- sigma
    vg
  },

  # ::modelo Legendre::
  # Pensado para distancias sobre la esfera: h se toma en kilómetros y se pasa
  # a un ángulo con el radio de la Tierra (6378.137 km), igual que en gstat
  legendre = function(h, rango = 0.5, sigma = 2, ...){
    stopifnot(
      'Model "legendre": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "legendre": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "legendre": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    angle <- h / (pi * 6378.137)
    vg <- sigma * (2 - (1 - rango^2) / (1 - 2 * rango * cos(angle) + rango^2))
    vg[h == 0] <- 0    # El variograma vale cero en el origen
    vg
  },

  # En teoría faltarían el modelo de intercepto y de error, pero por el
  # momento los voy a obviar: Los modelos de gstat y su implementación
  # en lo que sospecho es c se encuentran descritos aquí:
  # https://github.com/r-spatial/gstat/blob/main/src/vario_fn.c

  # -------------------------------------------- NON GSTAT MODELS ----------------------------------
  # ::modelo lanteujoul::
  lanteujoul = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "lanteujoul": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "lanteujoul": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "lanteujoul": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    sigma - sigma * sin(pi/2 * exp(- abs(h)/rango))
  },

  # ::modelo radon de orden 2::
  `radon ord 2` = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "radon ord 2": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "radon ord 2": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "radon ord 2": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    sigma-sigma*(1 + h/rango)*exp(-h/rango)
  },

  # ::modelo radon de orden 4::
  `radon ord 4` = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "radon ord 4": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "radon ord 4": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "radon ord 4": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    sigma-sigma*(1 + h/rango + (1/3) * (h/rango)^2) * exp(-h/rango)
  },

  # ::modelo cúbico:
  cubic = function(h, rango = 1, sigma = 2, ...){
    stopifnot(
      'Model "cubic": `h` must be numeric, non-negative and without missing values' =
        is.numeric(h) && !anyNA(h) && all(h >= 0),
      'Model "cubic": `rango` must be a single number greater than 0' =
        is.numeric(rango) && length(rango) == 1 && is.finite(rango) && rango > 0,
      'Model "cubic": `sigma` must be a single finite number' =
        is.numeric(sigma) && length(sigma) == 1 && is.finite(sigma)
    )
    vg     <- h
    bfr    <- vg < rango                  # antes de alcanzar el rango
    afr    <- !bfr                        # luego de alcanzar el rango

    vg[bfr] <- sigma*(7 * (h[bfr]/rango)**2 - 8.75 * (h[bfr]/rango)**3 + 3.5 * (h[bfr]/rango)**5 -
                        0.75 * (h[bfr]/rango) ** 7)
    vg[afr] <- sigma
    vg
  }
)

# FUNCIÓN PARA QUE EL USUARIO SEPA CUÁLES SON LAS POSIBLES OPCIONES DE MODELO

#' List or plot the available variogram models
#'
#' Shows the variogram models in [models], either as a grid of plots or as a
#' table with the arguments each model takes.
#'
#' @param txt If `FALSE` (the default), every model is plotted with
#'   `rango = 1.5`, `sigma = 5` and `alpha1 = 1.3`; dotted lines mark the range
#'   and the sill. If `TRUE`, the arguments of each model are printed instead.
#'
#' @return Called for its side effect. With `txt = TRUE`, the printed
#'   character matrix is returned invisibly.
#'
#' @examples
#' availablemodels(txt = TRUE)
#' availablemodels()
#'
#' @importFrom graphics abline curve par
#' @export
availablemodels <- function(txt = F){
  if (!txt){
    old = par(mfrow = c(5,5), mar = c(0.4, 0.4, 1.4, 0.4))
    invisible(
      lapply(names(models), \(model) {
        curve(models[[model]](x, rango = 1.5, sigma = 5, alpha1 = 1.3), 0, 8, ylab = "", xlab = "h",
              col = "darkblue",add = F, bty = 'n', ylim = c(0, 10),
              main = stringr::str_to_upper(model), lwd = 2, xaxt = 'n', yaxt = 'n')
        abline(v = 1.5, col = 'darkred', lty = 'dotted')
        abline(h = 5, col = 'darkred', lty = 'dotted')
      }))
    par(old)
  } else {
    arguments <- lapply(models, (\(x) paste0(names(formals(x)), collapse = ', '))) |>
      do.call(what = rbind)
    colnames(arguments) <- c('arguments')
    print(arguments)
  }
}
