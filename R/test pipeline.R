# ==================================================================================================
#                CASO DE USO - NUEVO FLUJO DE TRABAJO EN SpatFD (con traductor + compilación)
# ==================================================================================================

library(gstat)
library(SpatFD)
data(COKMexico)

# --------------------------------------------------------------------------------------------------
#                                PASO 0: Suavizamiento funcional y AFCP
# --------------------------------------------------------------------------------------------------
pm10 <- spatfd(Mex_PM10, coords = as.matrix(coord_PM10),
               basis  = "fourier", nbasis = 21, lambda = 0.000001,
               nharm  = 2, var    = 0.95,
               name   = "PM10")

no2 <- spatfd(NO2, coords = as.matrix(coord_NO2),
              basis  = "fourier", nbasis = 21, lambda = 0.000001,
              nharm  = 2, var    = 0.99,
              name   = "NO2")

pm10.no2 <- pm10 + no2

# --------------------------------------------------------------------------------------------------
#                                   PASO 1: Salida a gstat
# --------------------------------------------------------------------------------------------------
gstatPM10 <- togstat(pm10.no2, vars = 'PM10')
gstatALL  <- togstat(pm10.no2)

vgPM10 <- variogram(gstatPM10, cross = F)
vgALL  <- variogram(gstatALL)

# --------------------------------------------------------------------------------------------------
#                          MODELAMIENTO POR FUERA DEL PAQUETE (vía gstat)
# --------------------------------------------------------------------------------------------------

# --- Caso multivariado (crossed) ---
model <- vgm(647677.1, "Wav", 23317.05)
model <- vgm(127633,   "Wav",  9408.63, add.to = model)
model <- vgm(100, model = 'Nug', 0, add.to = model)

fitCrossed <- gstat::fit.lmc(v = vgALL, g = gstatALL,
                             model = model, fit.method = 6, correct.diagonal = 1.01)

# --- Caso univariado (lista de modelos, uno por PC) ---
model.PM10 <- list(
  'PM10.PC1' = vgm(psill = 20000, "Exp", range = 2103.25, nugget = 0),
  'PM10.PC2' = vgm(psill = 101494.96, "Exp", range = 1484.57, nugget = 0)
)

fitSingle <- lapply(names(gstatPM10$data), function(i) {
  gstat::fit.variogram(vgPM10[vgPM10$id == i, ], model = model.PM10[[i]], fit.method = 6)
})
names(fitSingle) <- names(model.PM10)

# --------------------------------------------------------------------------------------------------
#                     PASO 2: Traducir los modelos ajustados hacia el esquema SpatFD
# --------------------------------------------------------------------------------------------------
# toSpatFD() despacha según la clase del objeto de modelo entregado:
#   - fitCrossed  -> objeto 'gstat' con $model         -> toSpatFD.gstat()
#   - fitSingle   -> lista nombrada de 'variogramModel' -> toSpatFD.list()
# En ambos casos el resultado es el mismo esquema tabular: var, model, rango, sigma, alpha1

df_crossed <- toSpatFD(fitCrossed)
df_single  <- toSpatFD(fitSingle)

print(df_crossed)
print(df_single)

# --------------------------------------------------------------------------------------------------
#                          PASO 3: Compilar el df hacia funciones evaluables
# --------------------------------------------------------------------------------------------------
# compile_model() convierte cada variable del df en una única función de un argumento (h),
# suma de sus componentes, con los parámetros ya fijados.

compiled_crossed <- compile_model(df_crossed)
compiled_single  <- compile_model(df_single)

# --------------------------------------------------------------------------------------------------
#                          PASO 4: Verificación demostrativa
# --------------------------------------------------------------------------------------------------
# (a) Las funciones existen y son evaluables para cada variable/score:
str(compiled_crossed, max.level = 1)
str(compiled_single,  max.level = 1)

# (b) Evaluación numérica directa - confirma que no hay NA/Inf en un rango de distancias:
h_prueba <- seq(0, 30000, length.out = 50)

resultado_pm10_pc1 <- compiled_single[["PM10.PC1"]](h_prueba)
stopifnot(all(is.finite(resultado_pm10_pc1)))

# (c) Comparación visual contra el ajuste original de gstat, como verificación cruzada de que
#     el traductor + compile_model() reproduce el mismo modelo que gstat ajustó:

plot(h_prueba, resultado_pm10_pc1, type = 'l', col = 'darkblue', lwd = 2,
     main = "compile_model() vía toSpatFD() (SpatFD)",
     xlab = 'h', ylab = 'gamma',
     ylim = c(1, 1000000), xlim = c(0,20000))

plot(vgPM10[vgPM10$id == 'PM10.PC1', ], model = fitSingle$PM10.PC1,
     main = "Resultado gstat")
