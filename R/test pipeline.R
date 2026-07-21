# =============================================================================
# EJEMPLO DE USO - Nuevo pipeline de SpatFD
# =============================================================================
#
# El pipeline está diseñado con dos objetivos en mente:
#
#   1. SIMPLICIDAD: El flujo estándar (suavizamiento → variograma → ajuste →
#      kriging) se ejecuta con una cadena mínima de llamadas.
#
#   2. FLEXIBILIDAD: Cada paso produce un objeto intermedio accesible, por lo
#      que el usuario puede salirse del pipeline en cualquier momento para
#      aplicar otros métodos de modelamiento estadístico (e.g. modelos
#      no estacionarios, machine learning espacial, etc.).
#
# -----------------------------------------------------------------------------
# NOTA SOBRE EL CASO MULTIVARIADO
# -----------------------------------------------------------------------------
# Cuando se trabaja con varias variables, la forma recomendada es construir
# UN ÚNICO objeto SpatFD que las contenga a todas. Esto permite que
# varioSpatFD() calcule automáticamente los variogramas cruzados y que
# fitSpatFD() ajuste un Linear Model of Coregionalization (LMC) coherente.
#
# La forma de construirlo es con el operador +:
#
#   pm10.no2 <- pm10 + no2
#
# Aunque también es válido trabajar con cada variable por separado si solo
# se necesita el caso univariado (ver más abajo).
# =============================================================================

library(gstat)
library(SpatFD)
data(COKMexico)

# =============================================================================
# PASO 0: Suavizamiento funcional y AFCP
# =============================================================================
# spatfd() encapsula el suavizamiento (fda::smooth.basis) y el análisis
# funcional de componentes principales (fda::pca.fd) en un único objeto.
#
# El argumento `var = 0.95` activa la selección automática del número de
# componentes principales: se elige el mínimo de PCs que expliquen al menos
# el 95% de la varianza. Esto evita tener que determinar `nharm` a mano
# (nharm actúa aquí como tope máximo).

pm10 <- spatfd(Mex_PM10,
               coords = as.matrix(coord_PM10),
               basis  = "fourier",
               nbasis = 21,
               lambda = 0.000001,
               nharm  = 2,
               var    = 0.95,
               name   = "PM10")

no2 <- spatfd(NO2,
              coords = as.matrix(coord_NO2),
              basis  = "fourier",
              nbasis = 21,
              lambda = 0.000001,
              nharm  = 2,
              var    = 0.95,
              name   = "NO2")

# El operador + combina ambas variables en un único objeto multivariado.
# Las estaciones no comunes entre pm10 y no2 se incorporan automáticamente.
pm10.no2 <- pm10 + no2

pm10vg <- togstat(pm10.no2, vars='PM10')
no2vg <- togstat(pm10.no2, vars='NO2')
allvg <- togstat(pm10.no2)

allvg |> variogram() |> plot()
