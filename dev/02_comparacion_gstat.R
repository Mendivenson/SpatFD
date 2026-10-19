# ==================================================================================================
#                           SpatFD v2 - COMPARACIÓN DE RESULTADOS CON gstat
# ==================================================================================================
# La v2 ya no usa gstat para evaluar los variogramas ni para predecir: lo hace con código propio.
# Este script comprueba que los resultados son los mismos que da gstat:
#
#   A. Variogramas: el modelo guardado en el objeto es el mismo que ajustó gstat
#   B. Kriging:     scores predichos y varianzas, modelo sin variogramas cruzados
#   C. Cokriging:   scores predichos y varianzas, modelo lineal de corregionalización
#   D. Lambda:      curvas predichas y varianzas del método lambda univariado
#
# B y C se comparan con los dos tipos de kriging de gstat:
#   - simple   -> media conocida e igual a cero (beta = 0), lo que usa la versión CRAN
#   - ordinario -> media constante desconocida, lo que hace gstat por defecto; es el defecto de
#                  predict() en SpatFD
# D no tiene referencia directa en gstat: se compara con las fórmulas escritas a mano, evaluando
# las covarianzas con gstat. Ya no sigue a la versión CRAN; el cambio se informa aparte.

# Se parte de los objetos del flujo de uso: pm10.no2, los ajustes de gstat y las predicciones
source("dev/01_uso_paquete.R")

# Lo que se compara con gstat (variogramas compilados, scores predichos, curvas funcionales) no
# se le entrega al usuario: se saca de las funciones internas del paquete. Con load_all() están
# disponibles directamente; con el paquete instalado serían SpatFD:::get_models(), etc.

# Mayor diferencia entre dos resultados, relativa al tamaño de los valores de gstat
diferencia <- function(spatfd, gstat) max(abs(spatfd - gstat)) / max(abs(gstat))


# --------------------------------------------------------------------------------------------------
#                                       A. VARIOGRAMAS
# --------------------------------------------------------------------------------------------------
# Cada variograma guardado por add_model() se evalúa en las mismas distancias que el de gstat.

distancias <- seq(0, 30000, length.out = 200)

comparar_variogramas <- function(modelo, ajuste_gstat) {
  propios <- get_models(pm10.no2)[[modelo]]$variograms

  sapply(names(ajuste_gstat), function(id) {
    de_gstat <- variogramLine(ajuste_gstat[[id]], dist_vector = distancias)$gamma
    diferencia(propios[[id]](distancias), de_gstat)
  })
}

dif_variogramas <- list(
  "model 1"  = comparar_variogramas("model 1",  ajuste_pm10),
  "lmc wave" = comparar_variogramas("lmc wave", ajuste_lmc$model)
)

# Ejemplo gráfico: variograma del primer score de PM10
vario_pc1 <- vario_pm10[vario_pm10$id == "PM10.PC1", ]

plot(vario_pc1$dist, vario_pc1$gamma, pch = 19, col = "grey50",
     xlim = c(0, max(vario_pc1$dist)), ylim = c(0, max(vario_pc1$gamma) * 1.1),
     xlab = "Distancia", ylab = "Semivarianza", main = "Variograma de PM10.PC1")
lines(variogramLine(ajuste_pm10$PM10.PC1, dist_vector = distancias), col = "darkred", lwd = 5)
lines(distancias, get_models(pm10.no2)[["model 1"]]$variograms$PM10.PC1(distancias),
      col = "white", lwd = 1.5, lty = 2)
legend("bottomright", bty = "n",
       legend = c("Variograma empírico", "Modelo en gstat", "Modelo en SpatFD (línea blanca)"),
       pch = c(19, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 5, 1.5),
       col = c("grey50", "darkred", "grey30"))


# --------------------------------------------------------------------------------------------------
#                                         B. KRIGING
# --------------------------------------------------------------------------------------------------
# Referencia: gstat::krige() sobre cada score de PM10, con el modelo exponencial ajustado. Con
# beta = 0 es kriging simple; sin beta, ordinario.
#
# predict() solo entrega al usuario las curvas y su varianza. Los scores predichos, que es lo que
# se puede comparar con gstat, salen de la función interna krige_scores().

sitios_gstat <- sitios_nuevos
sp::coordinates(sitios_gstat) <- ~ X + Y

comparar_kriging <- function(tipo) {
  propio <- krige_scores(pm10.no2, sitios_nuevos, model = "model 1", kriging = tipo)

  sapply(names(ajuste_pm10), function(score) {
    beta <- if (tipo == "simple") 0 else numeric(0)
    de_gstat <- krige(reformulate("1", score), scores_pm10$data[[score]]$data, sitios_gstat,
                      model = ajuste_pm10[[score]], beta = beta, debug.level = 0)

    c(prediccion = diferencia(propio$pred[, score],     de_gstat$var1.pred),
      varianza   = diferencia(propio$variance[, score], de_gstat$var1.var))
  })
}

dif_kriging <- list(simple = comparar_kriging("simple"), ordinario = comparar_kriging("ordinary"))


# --------------------------------------------------------------------------------------------------
#                                        C. COKRIGING
# --------------------------------------------------------------------------------------------------
# Referencia: predict() de gstat sobre el modelo lineal de corregionalización. Para el cokriging
# simple se fija en cero la media de cada score; sin eso, gstat hace cokriging ordinario.

lmc_media_cero <- ajuste_lmc
for (score in names(lmc_media_cero$data)) lmc_media_cero$data[[score]]$beta <- 0

comparar_cokriging <- function(tipo) {
  modelo   <- if (tipo == "simple") lmc_media_cero else ajuste_lmc
  de_gstat <- as.data.frame(predict(modelo, sitios_gstat, debug.level = 0))
  propio   <- krige_scores(pm10.no2, sitios_nuevos, model = "lmc wave", kriging = tipo)

  sapply(colnames(propio$pred), function(score) {
    c(prediccion = diferencia(propio$pred[, score],     de_gstat[[paste0(score, ".pred")]]),
      varianza   = diferencia(propio$variance[, score], de_gstat[[paste0(score, ".var")]]))
  })
}

dif_cokriging <- list(simple = comparar_cokriging("simple"), ordinario = comparar_cokriging("ordinary"))


# --------------------------------------------------------------------------------------------------
#                                      D. MÉTODO LAMBDA
# --------------------------------------------------------------------------------------------------
# Referencia: las fórmulas escritas a mano, evaluando las covarianzas con gstat. Omega y c suman
# sobre los scores las covarianzas entre estaciones y con los sitios nuevos.
#
#   simple:    Omega pesos = c                    curva = media + pesos x (curvas - media)
#                                                 varianza = suma de mesetas - c' pesos
#   ordinario: | Omega 1 | |pesos|   |c|          curva = pesos x curvas (sin centrar: los
#              | 1'    0 | | mu  | = |1|          pesos suman 1 y la media se cancela)
#                                                 varianza = suma de mesetas - c' pesos - mu
#
# El cokriging funcional (lambda con variogramas cruzados entre variables) no tiene referencia en
# gstat: se verifica en tests/ contra el sistema por bloques escrito directamente.

estaciones <- pm10.no2$scores[!is.na(pm10.no2$scores[, "PM10.PC1"]), c("X", "Y")]
distancias_todas <- as.matrix(dist(rbind(estaciones, as.matrix(sitios_nuevos))))
entre_estaciones <- distancias_todas[rownames(estaciones), rownames(estaciones)]
a_sitios_nuevos  <- distancias_todas[rownames(estaciones), rownames(sitios_nuevos)]

covarianza <- function(modelo, d) variogramLine(modelo, dist_vector = d, covariance = TRUE)
omegas <- lapply(ajuste_pm10, covarianza, d = entre_estaciones)
cs     <- lapply(ajuste_pm10, covarianza, d = a_sitios_nuevos)
omega  <- Reduce(`+`, omegas)
c0     <- Reduce(`+`, cs)

tiempos   <- seq_len(nrow(pm10_mexico))
media     <- as.vector(fda::eval.fd(tiempos, get_pca(pm10.no2, "PM10")$meanfd))
observado <- fda::eval.fd(tiempos, get_fd(pm10.no2, "PM10"))
colnames(observado) <- colnames(get_fd(pm10.no2, "PM10")$coefs)
observado <- observado[, rownames(estaciones)]
mesetas   <- sum(sapply(ajuste_pm10, function(m) sum(m$psill)))

# Simple
pesos_simple  <- solve(omega, c0)
lambda_simple <- list(curvas   = media + (observado - media) %*% pesos_simple,
                      varianza = mesetas - colSums(c0 * pesos_simple))

# Ordinario
n <- nrow(estaciones)
sistema   <- rbind(cbind(omega, 1), c(rep(1, n), 0))
solucion  <- solve(sistema, rbind(c0, 1))
pesos_ord <- solucion[seq_len(n), ]
lambda_ordinario <- list(curvas   = observado %*% pesos_ord,
                         varianza = mesetas - colSums(c0 * pesos_ord) - solucion[n + 1, ])

pred_lambda_simple <- predict(pm10.no2, sitios_nuevos, model = "model 1", method = "lambda",
                              kriging = "simple")

dif_lambda <- list(
  simple    = c(curvas   = diferencia(pred_lambda_simple$PM10$pred, lambda_simple$curvas),
                varianza = diferencia(pred_lambda_simple$PM10$var,  lambda_simple$varianza)),
  ordinario = c(curvas   = diferencia(pred_lambda$PM10$pred, lambda_ordinario$curvas),
                varianza = diferencia(pred_lambda$PM10$var,  lambda_ordinario$varianza))
)

# REFERENCIA HISTÓRICA: la versión CRAN usaba kriging simple, aplicaba los pesos a las curvas sin
# centrar y tomaba la varianza de cada score de su valor propio. Solo se informa cuánto cambia
# frente al lambda simple de hoy; no debe coincidir.
valores_propios <- get_pca(pm10.no2, "PM10")$values
lambda_cran <- list(
  curvas   = observado %*% pesos_simple,
  varianza = Reduce(`+`, lapply(1:2, function(k) {
    valores_propios[k] - 2 * colSums(cs[[k]] * pesos_simple) +
      colSums(pesos_simple * (omegas[[k]] %*% pesos_simple))
  }))
)
cambio_cran <- c(curvas   = diferencia(pred_lambda_simple$PM10$pred, lambda_cran$curvas),
                 varianza = diferencia(pred_lambda_simple$PM10$var,  lambda_cran$varianza))


# --------------------------------------------------------------------------------------------------
#                                         RESULTADO
# --------------------------------------------------------------------------------------------------
# Diferencias relativas entre SpatFD y la referencia. Valores del orden de 1e-10 o menores son
# solo error numérico de redondeo: los resultados coinciden.

cat("\nA. Variogramas del modelo 'model 1'\n")
print(signif(dif_variogramas[["model 1"]], 2))

cat("\nA. Variogramas del modelo 'lmc wave' (directos y cruzados)\n")
print(signif(dif_variogramas[["lmc wave"]], 2))

for (tipo in c("simple", "ordinario")) {
  cat("\nB. Kriging", tipo, "de PM10\n")
  print(signif(dif_kriging[[tipo]], 2))
}

for (tipo in c("simple", "ordinario")) {
  cat("\nC. Cokriging", tipo, "de PM10 y NO2\n")
  print(signif(dif_cokriging[[tipo]], 2))
}

for (tipo in c("simple", "ordinario")) {
  cat("\nD. Método lambda", tipo, "de PM10\n")
  print(signif(dif_lambda[[tipo]], 2))
}

cat("\n   Lambda simple frente a la versión CRAN (curvas sin centrar, valores propios), solo informativo:\n")
print(signif(cambio_cran, 2))

tolerancia <- 1e-8
stopifnot(
  "Los variogramas no coinciden con gstat"    = all(unlist(dif_variogramas) < tolerancia),
  "El kriging no coincide con gstat"          = all(unlist(dif_kriging)     < tolerancia),
  "El cokriging no coincide con gstat"        = all(unlist(dif_cokriging)   < tolerancia),
  "El método lambda no coincide con la fórmula" = all(unlist(dif_lambda)    < tolerancia)
)
cat("\nTodo coincide (diferencia relativa menor que", tolerancia, ")\n")
