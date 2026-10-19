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
#   7. Predicción sobre un mapa
#   8. Comparación de modelos y métodos por validación cruzada
#
# La comprobación de que los resultados coinciden con gstat está en dev/02_comparacion_gstat.R

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

library(gstat)
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
#   - Como cualquier función de modelado de R, devuelve el objeto actualizado: hay que reasignarlo.
#   - Un objeto puede guardar varios modelos, para compararlos después.
#   - Si no se da un nombre, se nombran en secuencia: model 1, model 2, ...

pm10.no2 <- pm10.no2 |>
  add_model(ajuste_pm10) |>                          # queda como 'model 1'
  add_model(ajuste_lmc, name = "lmc wave")

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
#
# Hay dos métodos de predicción (argumento 'method'):
#   - "scores" (por defecto) -> predice los scores y reconstruye la curva con ellos
#   - "lambda"               -> la curva es la media más una combinación lineal de las curvas
#                               observadas (centradas). Con variogramas cruzados entre variables
#                               es cokriging funcional: las curvas de NO2 ayudan a predecir PM10.
#
# Y dos tipos de kriging (argumento 'kriging'):
#   - "ordinary" (por defecto) -> la media es una constante desconocida. Los pesos de cada variable
#                                 suman 1, así que da lo mismo centrar o no las curvas, y se pueden
#                                 usar modelos sin meseta (linear, power, ...).
#   - "simple"                 -> la media es conocida e igual a cero, como la de los scores.
#
# Para combinar las curvas de dos variables, lambda necesita que compartan los tiempos. Aquí PM10
# tiene 4344 instantes y NO2 4292, ambos con t = 1, 2, ...: lambda supone que empiezan a la misma
# hora y predice sobre el intervalo común [1, 4292], con un mensaje. Si los tiempos reales fueran
# otros, habría que darlos en 't' al crear los objetos con spatfd().

sitios_nuevos <- expand.grid(X = seq(478000, 505000, length.out = 3),
                             Y = seq(2135000, 2175000, length.out = 3))
rownames(sitios_nuevos) <- paste0("S", seq_len(nrow(sitios_nuevos)))

pred_kriging   <- predict(pm10.no2, sitios_nuevos, model = "model 1")
pred_cokriging <- predict(pm10.no2, sitios_nuevos, model = "lmc wave")
pred_lambda    <- predict(pm10.no2, sitios_nuevos, model = "model 1", method = "lambda")
pred_cok_lamb  <- predict(pm10.no2, sitios_nuevos, model = "lmc wave", method = "lambda")
pred_simple    <- predict(pm10.no2, sitios_nuevos, model = "model 1", kriging = "simple")

# Si el modelo no sirve para lo que se pide, el error muestra la tabla de modelos:
try(predict(pm10.no2, sitios_nuevos))                                        # falta escoger uno
try(predict(pm10.no2, sitios_nuevos, model = "model 1", variable = "NO2"))   # no cubre NO2


# --------------------------------------------------------------------------------------------------
#                                        6. RESULTADOS
# --------------------------------------------------------------------------------------------------
# Toda predicción es una lista con un elemento por variable:  <predicción>$<variable>
#   $pred  curvas predichas: una matriz con una fila por tiempo y una columna por sitio nuevo,
#          igual que los datos de entrada (no hace falta saber nada del paquete fda)
#   $var   varianza de predicción de la curva en cada sitio nuevo

pred_kriging                  # Resumen: sitios, variables y varianza de predicción
names(pred_cokriging)         # El cokriging predice PM10 y NO2

pred_cokriging$NO2$pred[1:5, ]   # Curvas predichas de NO2 (primeros 5 tiempos)
pred_cokriging$NO2$var           # Su varianza de predicción

# plot() dibuja las curvas predichas (en color) sobre las curvas observadas (en gris)
par(mfrow = c(1, 2))
plot(pred_kriging,   variable = "PM10", lwd = 0.5)                     # Kriging
plot(pred_cokriging, variable = "PM10", lwd = 0.5, col = "darkred")    # Cokriging
par(mfrow = c(1, 1))

# Los dos métodos con el mismo modelo: scores (kriging) y lambda
par(mfrow = c(1, 2))
plot(pred_kriging, variable = "PM10", lwd = 0.5, main = "PM10 - Método scores")
plot(pred_lambda,  variable = "PM10", lwd = 0.5, main = "PM10 - Método lambda", col = "darkgreen")
par(mfrow = c(1, 1))

# Cokriging funcional (lambda con el modelo lineal de corregionalización): solo cubre los tiempos
# comunes a PM10 y NO2
dim(pred_cok_lamb$PM10$pred)
plot(pred_cok_lamb, variable = "PM10", lwd = 0.5, main = "PM10 - Cokriging funcional", col = "purple")

# Comparación de modelos y métodos: varianza de predicción de PM10 en cada sitio nuevo
cbind(kriging           = pred_kriging$PM10$var,
      kriging_simple    = pred_simple$PM10$var,
      cokriging         = pred_cokriging$PM10$var,
      lambda            = pred_lambda$PM10$var,
      cokriging_lambda  = pred_cok_lamb$PM10$var)

# Con kriging ordinario se pueden usar modelos sin meseta, por ejemplo uno lineal con pepita
pm10.no2 <- add_model(pm10.no2, list(PM10.PC1 = vgm(30, "Lin", 0, nugget = 20000),
                                     PM10.PC2 = vgm(10, "Lin", 0, nugget = 5000)),
                      name = "lineal")
predict(pm10.no2, sitios_nuevos, model = "lineal")
try(predict(pm10.no2, sitios_nuevos, model = "lineal", kriging = "simple"))   # error: sin meseta


# --------------------------------------------------------------------------------------------------
#                                    7. PREDICCIÓN SOBRE UN MAPA
# --------------------------------------------------------------------------------------------------
# Si en lugar de sitios se entrega un mapa (objeto sf), predict() predice en una grilla regular
# de puntos dentro del mapa. plot() dibuja dos mapas lado a lado: la predicción y su varianza.
#   - n     -> número aproximado de puntos de la grilla
#   - Los puntos blancos son las estaciones donde se midió la variable.

pred_mapa <- predict(pm10.no2, map_mexico, model = "lmc wave", variable = "PM10", n = 15000)
plot(pred_mapa)

# Una curva solo se puede mostrar en un mapa en un instante: se escoge con 'time'.
# La varianza es la de toda la curva, así que su mapa no cambia.
plot(pred_mapa, time = 1000, main = "PM10")

# Las mismas predicciones se pueden ver como curvas
plot(pred_mapa, mode = "curves", lwd = 0.2)


# --------------------------------------------------------------------------------------------------
#                    8. COMPARACIÓN DE MODELOS Y MÉTODOS POR VALIDACIÓN CRUZADA
# --------------------------------------------------------------------------------------------------
# leave_one_out() deja por fuera una estación a la vez, repite el análisis sin ella (suavizamiento
# y componentes principales), predice su curva con las demás y mide qué tan lejos queda de la
# curva observada. Los modelos de variograma no se vuelven a ajustar.
#
# Por defecto compara cada modelo que cubre la variable con cada método que admite.

validacion <- leave_one_out(pm10.no2, variable = "PM10")
validacion                    # Tabla de comparación, del menor al mayor error medio

validacion$errors             # Error de cada estación con cada combinación

# También se pueden escoger las combinaciones con una tabla de modelo y método
leave_one_out(pm10.no2, variable = "PM10",
              design = data.frame(model  = c("model 1", "lmc wave"),
                                  method = c("scores",  "scores")))

# Curvas predichas al dejar por fuera cada estación (en color) frente a las observadas (en gris).
# Son matrices tiempos x estaciones, así que se grafican con matplot().
matplot(pm10_mexico, type = "l", lty = 1, col = "grey70",
        xlab = "Tiempo", ylab = "PM10", main = "PM10 - Validación cruzada (lmc wave)")
matlines(validacion$pred[["lmc wave / scores"]], lty = 1, lwd = 0.5, col = "darkred")
