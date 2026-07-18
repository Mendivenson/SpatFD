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

# =============================================================================
# PASO 1: Variograma empírico
# =============================================================================
# varioSpatFD() detecta automáticamente el caso univariado o multivariado:
#   - Si el objeto tiene una sola variable → variograma directo.
#   - Si tiene varias → variogramas directos + cruzados (covariograma cruzado).
#
# El argumento `vars` permite calcular el variograma solo para un subconjunto
# de variables. Esto es útil para explorar la estructura espacial de cada
# variable antes de ajustar el LMC conjunto.

# Caso univariado: solo PM10 (a partir del objeto univariado original)
varioSpatFD(pm10) |> plot(main = 'PM10')

# Exploración univariada sobre el objeto multivariado:
# útil para afinar el modelo de PM10 antes del ajuste conjunto
varioSpatFD(pm10.no2, vars = c('PM10', 'NO2')) |> plot(main = 'PM10 x NO2')

# Objetos intermedios: el campo $variogram expone el data.frame estándar de
# gstat, compatible con cualquier función del ecosistema de R.
vsfd_pm10 <- varioSpatFD(pm10.no2, vars = "PM10")   # solo variograma de PM10
vsfd_all  <- varioSpatFD(pm10.no2)                  # variogramas + cruzados

# Ejemplo de acceso directo al variograma empírico para uso externo:
head(vsfd_pm10$variogram)  # data.frame con columnas np, dist, gamma, id, ...

# =============================================================================
# PASO 2: Ajuste del modelo teórico de variograma
# =============================================================================
# fitSpatFD() usa la misma llamada para el caso univariado y multivariado.
# La detección del caso es interna.
#
# CASO MULTIVARIADO (LMC): se pasa un único variogramModel como estructura
# base. fitSpatFD() llama a gstat::fit.lmc() internamente.

model_lmc <- vgm(647677.1, "Gau", 23317.05)
model_lmc <- vgm(127633,   "Wav",  9408.63, add.to = model_lmc)

fitsfd_all <- fitSpatFD(v = vsfd_all, model = model_lmc)

# CASO UNIVARIADO: se puede pasar una lista con un modelo por PC.
# Si se pasa menos modelos que PCs, fill = TRUE los recicla automáticamente.
modelos_pm10 <- list(
  vgm(psill = 2634000,   "Exp", range =  2103.25, nugget = 0),
  vgm(psill =  101495,   "Exp", range =  1484.57, nugget = 0),
  vgm(psill =   53673,   "Exp", range = 42406.00, nugget = 0)
)

fitsfd_pm10 <- fitSpatFD(v = vsfd_pm10, model = modelos_pm10)

# Los objetos fit también son accesibles directamente:
#   fitsfd_pm10$fit  → data.frame con el modelo ajustado por PC
#   fitsfd_all$fit   → objeto gstat con el LMC ajustado
# Esto permite, por ejemplo, extraer los parámetros para un informe o
# usar el objeto gstat directamente con predict.gstat().

fitsfd_pm10 |> plot()
fitsfd_all  |> plot()

# La predicción simplemente se hace prediciendo los scores y reutilizando la función 
# recons_fd

p_all <- predict(fitsfd_all, newcoords = matrix(c(300000, 700000), nrow = 1))

recons_fd(p_all, 'PM10')
recons_fd(p_all, "PM10")
recons_fd(p_all, c("PM10", "NO2"))

plot(recons_fd(p_all))
plot(recons_fd(p_all), together = FALSE)

fda::eval.fd(1:365, recons_fd(p_all, "PM10")[["PM10"]])
