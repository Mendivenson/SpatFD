models <- list(

  # ---> Efecto pepita
  nugget       = function(h, rango = 1, sigma = 2, ...){
    rep(sigma, length(h))
  },

  # ---> Modelo triangular
  triangular   = function(h, rango = 1, sigma = 2, ...){
    vg     <- h
    # Caso 1 (c1): Antes de alcanzar el rango, Caso 2 (c2): Después del rango
    c1     <- vg < rango
    c2     <- !c1
    # Cálculo de la covarianza
    vg[c1] <- sigma * h[c1] / rango
    vg[c2] <- sigma
    vg
  },

  # ---> Modelo circular
  circular    = function(h, rango = 1, sigma = 2, ...){
    vg     <- h
    # Caso 1 (c1): Antes de alcanzar el rango, Caso 2 (c2): Después del rango
    c1    <- vg < rango
    c2    <- !c1
    # Cálculo de la covarianza
    vg[c1] <- sigma * (1.5*(h[c1]/rango) - 0.5*((h[c1]**3)/(rango**3)))
    vg[c2] <- sigma
    vg
  },

  # ---> Modelo esférico
  spherical   = function(h, rango = 1, sigma = 2, ...){
    vg     <- h
    # Caso 1 (c1): Antes de alcanzar el rango, Caso 2 (c2): Después del rango
    c1    <- vg < rango
    c2    <- !c1
    # Cálculo de la covarianza
    vg[c1] <- sigma - sigma * (1 - 1.5 * h[c1]/rango + 0.5 * (h[c1]/rango)^3)
    vg[c2] <- sigma
    vg
  },

  # ---> Modelo cúbico
  cubical     = function(h, rango = 1, sigma = 2, ...){
    vg     <- h
    # Caso 1 (c1): Antes de alcanzar el rango, Caso 2 (c2): Después del rango
    c1    <- vg < rango
    c2    <- !c1
    # Cálculo de la covarianza
    vg[c1] <- sigma*(7 * (h[c1]/rango)**2 - 8.75 * (h[c1]/rango)**3 + 3.5 * (h[c1]/rango)**5 - 0.75 * (h[c1]/rango) ** 7)
    vg[c2] <- sigma
    vg
  },

  # ---> Modelo exponencial
  exponential = function(h, rango = 1, sigma = 2, ...){
    vg  <- h
    # Cálculo de la covarianza
    vg  <- sigma - sigma * exp(- h/rango)
    vg
  },
  wave        = function(h, rango = 1, sigma = 2, ...){
    vg  <- h
    # Cálculo de la covarianza
    vg  <- sigma-sigma*(rango/h)*sin(h/rango)
    vg
  },

  # ---> Modelo pentamodel
  pentamodel  = function(h, rango = 1, sigma = 2, ...){
    vg     <- h
    # Caso 1 (c1): Antes de alcanzar el rango, Caso 2 (c2): Después del rango
    c1     <- vg < rango
    c2     <- !c1
    # Cálculo de la covarianza
    vg[c1] <- sigma * ( (22/3) * (h[c1]/rango)^2 - 33 * (h[c1]/rango)^4 + 38.5 * (h[c1]/rango)^5 -
                        16.5 * (h[c1]/rango)^7 + 5.5 * (h[c1]/rango)^9 - (5/6) * (h[c1]/rango)^11)
    vg[c2] <- sigma
    vg
  },

  # ---> Modelo radon 2
  `radon ord 2`   = function(h, rango = 1, sigma = 2, ...){
    vg  <- h
    # Cálculo de la covarianza
    vg  <- sigma-sigma*(1 + h/rango)*exp(-h/rango)
    vg
  },

  # ---> Modelo radon 4
  `radon ord 4`   = function(h, rango = 1, sigma = 2, ...){
    vg  <- h
    # Cálculo de la covarianza
    vg  <- sigma-sigma*(1 + h/rango + (1/3) * (h/rango)^2) * exp(-h/rango)
    vg
  },

  # ---> Modelo gaussiano
  gaussian   = function(h, rango = 1, sigma = 2, ...){
    vg  <- h
    # Cálculo de la covarianza
    vg  <- sigma - sigma * exp(- (h/rango)^2)
    vg
  },

  # ---> Modelo Lanteujoul
  lanteujoul   = function(h, rango = 1, sigma = 2, ...){
    vg  <- h
    # Cálculo de la covarianza
    vg  <- sigma - sigma * sin(pi/2 * exp(- abs(h)/rango))
    vg
  },

  # --> Modelo matern
  matern      = function(h, rango = 1, sigma = 2, alpha1 = 0.25, ...){
    vg <- h
    # Cálculo de la covarianza
    vg <- sigma - sigma *(1 / ((2**(alpha1-1)) * gamma(alpha1)))*((h/rango)**alpha1)*besselK(x = h/rango,nu = alpha1)
    vg
  }
)


class(models) <- 'models.SpatFD'

par(mfrow = c(4,4), mar = c(0.2, 0.2, 1.4, 0.2))
invisible(
  lapply(names(models), \(model) {
    curve(models[[model]](x, rango = 5, sigma = 5, alpha = 0.25), 0, 70, ylab = "", xlab = "h",
          col = "darkblue",add = F,
          ylim = c(0,7), main = stringr::str_to_upper(model), lwd = 2, xaxt = 'n', yaxt = 'n')
  }))
