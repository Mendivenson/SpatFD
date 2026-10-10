# ==================================================================================================
#                              SpatFD v2 - FLUJO DE USO DEL PAQUETE
# ==================================================================================================
# De las curvas observadas a las curvas predichas en sitios sin medición:
#
#   1. Datos
#   2. Objeto SpatFD: suavizamiento y componentes principales funcionales
#   3. Variogramas de los scores (por fuera del paquete, aquí con gstat)
#   4. Agregar los modelos ajustados al objeto
#   5. Predicción en sitios nuevos
#   6. Resultados
#
# La comprobación de que los resultados coinciden con gstat está en dev/02_comparacion_gstat.R

library(gstat)

# --------------------------------------------------------------------------------------------------
#                                          1. DATOS
# --------------------------------------------------------------------------------------------------
# Contaminación del aire en Ciudad de México: PM10 y NO2 medidos a lo largo del tiempo en
# varias estaciones. Los datos vienen incluidos en el paquete (ver ?air_mexico).
#   - pm10_mexico, no2_mexico                -> una columna por estación, una fila por momento de medición
#   - coords_pm10_mexico, coords_no2_mexico  -> coordenadas de las estaciones

devtools::load_all()
data(air_mexico)

# --------------------------------------------------------------------------------------------------
#                                     2. OBJETO SpatFD
# --------------------------------------------------------------------------------------------------
# spatfd() suaviza las curvas y calcula sus componentes principales funcionales. De aquí en
# adelante cada curva queda resumida por sus scores, que es lo que se modela espacialmente.

pm10 <- spatfd(pm10_mexico, coords = as.matrix(coords_pm10_mexico),
               basis = "fourier", nbasis = 21, lambda = 0.000001,
               nharm = 2, var = 0.95,
               name  = "PM10")

no2 <- spatfd(no2_mexico, coords = as.matrix(coords_no2_mexico),
              basis = "fourier", nbasis = 21, lambda = 0.000001,
              nharm = 2, var = 0.99,
              name  = "NO2")

# Las variables se combinan con '+'. Cada una conserva sus propias estaciones.
pm10.no2 <- pm10 + no2

summary(pm10.no2)                                  # Componentes y varianza explicada
plot(pm10.no2, variable = "PM10")                  # Curvas suavizadas
plot(pm10.no2, variable = "PM10", mode = "pca")    # Componentes principales


# --------------------------------------------------------------------------------------------------
#                              3. VARIOGRAMAS DE LOS SCORES (gstat)
# --------------------------------------------------------------------------------------------------
# El paquete no ajusta variogramas: entrega los scores y el usuario los modela con la herramienta
# que prefiera. togstat() los deja listos para gstat.

scores_pm10  <- togstat(pm10.no2, vars = "PM10")
scores_todos <- togstat(pm10.no2)

vario_pm10  <- variogram(scores_pm10, cross = FALSE)
vario_todos <- variogram(scores_todos)

# Modelo A: un variograma exponencial por cada score de PM10 (sin variogramas cruzados)
inicial_pm10 <- list(
  "PM10.PC1" = vgm(psill = 20000,     "Exp", range = 2103.25, nugget = 0),
  "PM10.PC2" = vgm(psill = 101494.96, "Exp", range = 1484.57, nugget = 0)
)

ajuste_pm10 <- lapply(names(inicial_pm10), function(score) {
  fit.variogram(vario_pm10[vario_pm10$id == score, ], model = inicial_pm10[[score]], fit.method = 6)
})
names(ajuste_pm10) <- names(inicial_pm10)

# Modelo B: modelo lineal de corregionalización para PM10 y NO2 (con variogramas cruzados)
inicial_lmc <- vgm(647677.1, "Wav", 23317.05)
inicial_lmc <- vgm(127633,   "Wav",  9408.63, add.to = inicial_lmc)
inicial_lmc <- vgm(100,      "Nug",        0, add.to = inicial_lmc)

ajuste_lmc <- fit.lmc(v = vario_todos, g = scores_todos,
                      model = inicial_lmc, fit.method = 6, correct.diagonal = 1.01)


# --------------------------------------------------------------------------------------------------
#                              4. AGREGAR LOS MODELOS AL OBJETO
# --------------------------------------------------------------------------------------------------
# add_model() guarda el modelo dentro del objeto SpatFD, listo para predecir.
#   - Recibe directamente lo que devuelve gstat.
#   - Modifica el objeto en el sitio: no hay que reasignarlo.
#   - Un objeto puede guardar varios modelos, para compararlos después.
#   - Si no se da un nombre, se nombran en secuencia: model 1, model 2, ...

add_model(pm10.no2, ajuste_pm10)                     # queda como 'model 1'
add_model(pm10.no2, ajuste_lmc, name = "lmc wave")

# Tabla resumen de los modelos del objeto
pm10.no2$models

# Opcional: toSpatFD() muestra cómo queda un modelo en el formato propio del paquete
toSpatFD(ajuste_pm10)


# --------------------------------------------------------------------------------------------------
#                                 5. PREDICCIÓN EN SITIOS NUEVOS
# --------------------------------------------------------------------------------------------------
# predict() recibe los sitios nuevos y el nombre del modelo. Predice los scores por kriging
# simple y reconstruye las curvas. El modelo decide qué se resuelve:
#   - Sin variogramas cruzados -> kriging, un score a la vez
#   - Con variogramas cruzados -> cokriging, todos los scores en conjunto

sitios_nuevos <- expand.grid(X = seq(478000, 505000, length.out = 3),
                             Y = seq(2135000, 2175000, length.out = 3))
rownames(sitios_nuevos) <- paste0("S", seq_len(nrow(sitios_nuevos)))

pred_kriging   <- predict(pm10.no2, sitios_nuevos, model = "model 1")
pred_cokriging <- predict(pm10.no2, sitios_nuevos, model = "lmc wave")

# Si el modelo no sirve para lo que se pide, el error muestra la tabla de modelos:
try(predict(pm10.no2, sitios_nuevos))                                        # falta escoger uno
try(predict(pm10.no2, sitios_nuevos, model = "model 1", variable = "NO2"))   # no cubre NO2


# --------------------------------------------------------------------------------------------------
#                                        6. RESULTADOS
# --------------------------------------------------------------------------------------------------
# Toda predicción es una lista con un elemento por variable:  <predicción>$<variable>
#   $pred  curvas predichas (objeto fd del paquete fda), una por sitio nuevo
#   $var   varianza de predicción de la curva en cada sitio nuevo

pred_kriging                  # Resumen: sitios, variables y varianza de predicción
names(pred_cokriging)         # El cokriging predice PM10 y NO2

pred_cokriging$NO2$pred       # Curvas predichas de NO2
pred_cokriging$NO2$var        # Su varianza de predicción

# Curvas predichas (en color) sobre las curvas observadas (en gris)
par(mfrow = c(1, 2))

plot(attr(pm10.no2, "fd")$PM10, col = "grey80", lty = 1, main = "PM10 - Kriging")
lines(pred_kriging$PM10$pred, col = "darkblue", lty = 1, lwd = 0.5)

plot(attr(pm10.no2, "fd")$PM10, col = "grey80", lty = 1, main = "PM10 - Cokriging")
lines(pred_cokriging$PM10$pred, col = "darkred", lty = 1, lwd = 0.5)

par(mfrow = c(1, 1))

# Comparación de los dos modelos: varianza de predicción de PM10 en cada sitio nuevo
cbind(kriging   = pred_kriging$PM10$var,
      cokriging = pred_cokriging$PM10$var)
