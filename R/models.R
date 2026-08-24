models <- list(

  # ----------------------------------- GSTAT SUPPORTED MODELS -------------------------------------

  # ::efecto pepita::
  nugget = function(h, sigma = 2, ...){
    rep(sigma, length(h))
  },

  # ::modelo exponencial::
  exponential = function(h, rango = 1, sigma = 2, ...){
    sigma - sigma * exp(- h/rango)
  },

  # ::modelo esférico::
  spherical = function(h, rango = 1, sigma = 2, ...){
    vg     <- h
    bfr    <- vg < rango                   # Antes de alcanzar el rango
    afr    <- !bfr                         # Después de alcanzar el rango

    vg[bfr] <- sigma - sigma * (1 - 1.5 * h[bfr]/rango + 0.5 * (h[bfr]/rango)^3)
    vg[afr] <- sigma
    vg
  },

  # ::modelo gaussiano::
  gaussian = function(h, rango = 1, sigma = 2, ...){
    sigma - sigma * exp(- ((h/rango)^2))
  },

  # ::modelo estable::
  stable = function(h, rango = 1, sigma = 2, alpha1 = 1.3,...){
    vg  <- h
    vg  <- sigma - sigma * exp( -((h/rango)^alpha1))
    vg
  },

  # ::modelo matern::
  matern = function(h, rango = 1, sigma = 2, alpha1 = 0.25, ...){
    vg <- sigma - sigma *
      (1 / ((2**(alpha1-1)) * gamma(alpha1))) * ((h/rango)**alpha1) *
      besselK(x = h/rango, nu = alpha1)
    vg[h == 0] <- 0         # EVITAR: Nan en Bessel eval
    vg
  },

  # ::modelo matern (parametrización Stein)::
  `stein's matern` = function(h, rango = 1, sigma = 2, alpha1 = 0.25, ...){
    arg <- 2 * sqrt(alpha1) * h / rango
    vg <- sigma - sigma *
      (2^(1 - alpha1) / gamma(alpha1)) * (arg^alpha1) *
      besselK(x = arg, nu = alpha1)
    vg[h == 0] <- 0        # EVITAR: Nan en Bessel eval
    vg
  },

  # ::modelo circular::
  circular = function(h, rango = 1, sigma = 2, ...){
    vg  <- h
    bfr <- vg < rango       # Antes de alcanzar el rango
    afr <- !bfr             # Después de alcanzar el rango

    vg[bfr] <- sigma * (1.5*(h[bfr]/rango) - 0.5*((h[bfr]**3)/(rango**3)))
    vg[afr] <- sigma
    vg
  },

  # ::modelo lineal::
  linear = function(h, sigma = 2, ...){
    sigma * h
  },

  # ::modelo Bessel::
  bessel = function(h, rango = 1, sigma = 2, alpha1 = 0.25, ...){
    vg <- sigma - sigma *
      (2**((alpha1/2)-1)) * gamma(alpha1/2) * ((h/rango)**(1-(alpha1/2))) *
      besselJ(x = h/rango, nu = ((alpha1/2)-1))
    vg[h == 0] <- 0
    vg
  },

  # ::modelo pentaesférico::
  pentamodel = function(h, rango = 1, sigma = 2, ...){
    vg  <- h
    bfr <- vg < rango          # Antes de alcanzar el rango
    afr <- !bfr                # Luego de alcanzar el rango

    vg[bfr] <- sigma * ( (22/3) * (h[bfr]/rango)^2 - 33 * (h[bfr]/rango)^4 + 38.5 * (h[bfr]/rango)^5 -
                          16.5 * (h[bfr]/rango)^7 + 5.5 * (h[bfr]/rango)^9 - (5/6) * (h[bfr]/rango)^11)
    vg[afr] <- sigma
    vg
  },

  # ::modelo periódico::
  periodic = function(h, rango = 1, sigma = 2, ...){
    sigma - sigma * cos((2 * pi * h) / rango)
  },

  # ::modelo wave::
  wave        = function(h, rango = 1, sigma = 2, ...){
    vg  <- h

    vg  <- sigma - sigma * (rango/ (pi * h)) * sin(pi * h/rango)
    vg[h == 0] <- 0
    vg
  },

  # ::modelo hueco::
  hole = function(h, rango = 1, sigma = 2, ...){
    vg <- sigma - sigma * sin(h/rango) / (h/rango)
    vg[h == 0] <- 0    # EVITAR: 0/0
    vg
  },

  # ::modelo logarítmico::
  logarithmic = function(h,sigma = 2, ...){
    vg <- sigma * log(h)
    vg[h == 0] <- 0    # EVITAR: -Inf log eval
    vg
  },

  # ::modelo power::
  power = function(h, rango = 1, sigma = 2, ...){
    sigma * (h ^ rango)
  },

  # ::modelo spline::
  spline = function(h, rango = 1, sigma = 2, ...){
    vg <- h * h * log(h)
    vg[h == 0] <- 0    # EVITAR: -Inf log eval
    vg
  },

  # ::modelo Legendre::
  legendre = function(h, rango = 1, sigma = 2, ...){
    vg <- h
    vg
  },

  # En teoría faltarían el modelo de intercepto y de error, pero por el
  # momento los voy a obviar: Los modelos de gstat y su implementación
  # en lo que sospecho es c se encuentran descritos aquí:
  # https://github.com/r-spatial/gstat/blob/main/src/vario_fn.c

  # -------------------------------------------- NON GSTAT MODELS ----------------------------------
  # ::modelo lanteujoul::
  lanteujoul = function(h, rango = 1, sigma = 2, ...){
    sigma - sigma * sin(pi/2 * exp(- abs(h)/rango))
  },

  # ::modelo radon de orden 2::
  `radon ord 2` = function(h, rango = 1, sigma = 2, ...){
    sigma-sigma*(1 + h/rango)*exp(-h/rango)
  },

  # ::modelo radon de orden 4::
  `radon ord 4` = function(h, rango = 1, sigma = 2, ...){
    sigma-sigma*(1 + h/rango + (1/3) * (h/rango)^2) * exp(-h/rango)
  },

  # ::modelo cúbico:
  cubic = function(h, rango = 1, sigma = 2, ...){
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
