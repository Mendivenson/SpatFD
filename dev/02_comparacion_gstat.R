# ==================================================================================================
#                           SpatFD v2 - COMPARACIÓN DE RESULTADOS CON gstat
# ==================================================================================================
# La v2 ya no usa gstat para evaluar los variogramas ni para predecir: lo hace con código propio.
# Este script comprueba que los resultados son los mismos que da gstat:
#
#   A. Variogramas: el modelo guardado en el objeto es el mismo que ajustó gstat
#   B. Kriging:     scores predichos y varianzas, modelo sin variogramas cruzados
#   C. Cokriging:   scores predichos y varianzas, modelo lineal de corregionalización
#
# La referencia es el kriging simple de gstat con media cero (beta = 0), que es lo que usa la
# versión CRAN de SpatFD por dentro.

# Se parte de los objetos del flujo de uso: pm10.no2, los ajustes de gstat y las predicciones
source("dev/01_uso_paquete.R")

# Mayor diferencia entre dos resultados, relativa al tamaño de los valores de gstat
diferencia <- function(spatfd, gstat) max(abs(spatfd - gstat)) / max(abs(gstat))


# --------------------------------------------------------------------------------------------------
#                                       A. VARIOGRAMAS
# --------------------------------------------------------------------------------------------------
# Cada variograma guardado por add_model() se evalúa en las mismas distancias que el de gstat.

distancias <- seq(0, 30000, length.out = 200)

comparar_variogramas <- function(modelo, ajuste_gstat) {
  propios <- attr(pm10.no2, "models")[[modelo]]$variograms

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
lines(distancias, attr(pm10.no2, "models")[["model 1"]]$variograms$PM10.PC1(distancias),
      col = "white", lwd = 1.5, lty = 2)
legend("bottomright", bty = "n",
       legend = c("Variograma empírico", "Modelo en gstat", "Modelo en SpatFD (línea blanca)"),
       pch = c(19, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 5, 1.5),
       col = c("grey50", "darkred", "grey30"))


# --------------------------------------------------------------------------------------------------
#                                         B. KRIGING
# --------------------------------------------------------------------------------------------------
# Referencia: gstat::krige() sobre cada score de PM10, con el modelo exponencial ajustado.

#
# predict() solo entrega al usuario las curvas y su varianza. Los scores predichos, que es lo que
# se puede comparar con gstat, salen de la función interna krige_scores().

sitios_gstat <- sitios_nuevos
sp::coordinates(sitios_gstat) <- ~ X + Y

kriging_spatfd <- krige_scores(pm10.no2, sitios_nuevos, model = "model 1")

dif_kriging <- sapply(names(ajuste_pm10), function(score) {
  kriging_gstat <- krige(reformulate("1", score), scores_pm10$data[[score]]$data, sitios_gstat,
                         model = ajuste_pm10[[score]], beta = 0, debug.level = 0)

  c(prediccion = diferencia(kriging_spatfd$pred[, score],     kriging_gstat$var1.pred),
    varianza   = diferencia(kriging_spatfd$variance[, score], kriging_gstat$var1.var))
})


# --------------------------------------------------------------------------------------------------
#                                        C. COKRIGING
# --------------------------------------------------------------------------------------------------
# Referencia: predict() de gstat sobre el modelo lineal de corregionalización. Para que sea
# cokriging simple se fija en cero la media de cada score.

lmc_media_cero <- ajuste_lmc
for (score in names(lmc_media_cero$data)) lmc_media_cero$data[[score]]$beta <- 0

cokriging_gstat <- as.data.frame(predict(lmc_media_cero, sitios_gstat, debug.level = 0))

cokriging_spatfd <- krige_scores(pm10.no2, sitios_nuevos, model = "lmc wave")

dif_cokriging <- sapply(colnames(cokriging_spatfd$pred), function(score) {
  c(prediccion = diferencia(cokriging_spatfd$pred[, score],     cokriging_gstat[[paste0(score, ".pred")]]),
    varianza   = diferencia(cokriging_spatfd$variance[, score], cokriging_gstat[[paste0(score, ".var")]]))
})


# --------------------------------------------------------------------------------------------------
#                                         RESULTADO
# --------------------------------------------------------------------------------------------------
# Diferencias relativas entre SpatFD y gstat. Valores del orden de 1e-10 o menores son solo
# error numérico de redondeo: los resultados coinciden.

cat("\nA. Variogramas del modelo 'model 1'\n")
print(signif(dif_variogramas[["model 1"]], 2))

cat("\nA. Variogramas del modelo 'lmc wave' (directos y cruzados)\n")
print(signif(dif_variogramas[["lmc wave"]], 2))

cat("\nB. Kriging de PM10\n")
print(signif(dif_kriging, 2))

cat("\nC. Cokriging de PM10 y NO2\n")
print(signif(dif_cokriging, 2))

tolerancia <- 1e-8
stopifnot(
  "Los variogramas no coinciden con gstat" = all(unlist(dif_variogramas) < tolerancia),
  "El kriging no coincide con gstat"       = all(dif_kriging     < tolerancia),
  "El cokriging no coincide con gstat"     = all(dif_cokriging   < tolerancia)
)
cat("\nTodo coincide con gstat (diferencia relativa menor que", tolerancia, ")\n")
