# ==================================================================================================
#                             CASO DE USO - NUEVO FLUJO DE TRABAJO EN SpatFD
# ==================================================================================================
#
# El pipeline se diseña con el objetivo en mente de permitir diversas herramientas de modelamiento
# (Fuera del paquete) proporcionando herramientas para sacar los datos de los scores hacia otros
# paquetes y para traer los modelos de vuelta hacia SpatFD o para definirlos de forma estándar de
# modo que aún cuando se usen métodos/modelos poco conocidos el análisis funcional aún se haga
# enteramentedentro del paquete.
#
# ==================================================================================================

library(gstat)
library(SpatFD)
data(COKMexico)

# --------------------------------------------------------------------------------------------------
#                                PASO 0: Suavizamiento funcional y AFCP
# --------------------------------------------------------------------------------------------------
# spatfd() encapsula el suavizamiento (fda::smooth.basis) y el análisis funcional de componentes
# principales (fda::pca.fd) en un único objeto.
#
# El argumento `var = 0.95` activa la selección automática del número de componentes principales:
# se elige el mínimo de PCs que expliquen al menos el 95% de la varianza. Esto evita tener que
# determinar `nharm` a mano (nharm actúa aquí como tope máximo).

pm10 <- spatfd(Mex_PM10, coords = as.matrix(coord_PM10),
               basis  = "fourier", nbasis = 21, lambda = 0.000001,
               nharm  = 2, var    = 0.95,
               name   = "PM10")

no2 <- spatfd(NO2, coords = as.matrix(coord_NO2),
              basis  = "fourier", nbasis = 21, lambda = 0.000001,
              nharm  = 2, var    = 0.99,
              name   = "NO2")

# El operador + combina ambas variables en un único objeto multivariado. Las estaciones no comunes
# entre pm10 y no2 se incorporan automáticamente.
pm10.no2 <- pm10 + no2

plot(pm10.no2, mode =  'functional', variable = 'PM10')
plot(pm10.no2, mode =  'functional', variable = 'NO2')

# --------------------------------------------------------------------------------------------------
#                                   PASO 1: Salida a otros paquetes
# --------------------------------------------------------------------------------------------------
# Por el momento, el paquete sólo ofrece modelamiento nativo en gstat, pero en el futuro se agregarán
# otros paquetes de modelamiento geoestadístico.
#
# La función gstat esencialmente combina en un solo objeto de gstat los diferentes scores con el fin
# de dejar listo el objeto para generar el variograma. En el caso univariado debe generarse el
# variograma con el argumento cross = F.
#
# NOTA: El flujo correcto de trabajo sería agregar todas las variables en un solo objeto y a partir
#       allí generar las combinaciones de variogramas necesarias para el modelamiento.

gstatPM10 <- togstat(pm10.no2, vars='PM10')
gstatNO2 <- togstat(pm10.no2, vars='NO2')
gstatALL <- togstat(pm10.no2)


vgPM10 <- variogram(gstatPM10, cross = F)
vgALL  <- variogram(gstatALL)

plot(vgPM10)
plot(vgALL)
# --------------------------------------------------------------------------------------------------
#                                MODELAMIENTO POR FUERA DEL PAQUETE
# --------------------------------------------------------------------------------------------------
# Idealmente se debe modelar por medio de modelos usuales de geoestadística. Sin embargo, también es
# posible modelar de forma no convencional y con modelos no convencionales, pero la inclusión de los
# modelos en el objeto SpatFD se vuelve un poco más compleja

# VÍA MODELOS CONVENCIONALES (gstat):
model <- vgm(647677.1, "Wav", 23317.05)
model <- vgm(127633,   "Wav",  9408.63, add.to = model)
model <- vgm(100, model = 'Nug', 0, add.to = model)

# Caso multivariado
fitCrossed <-  gstat::fit.lmc(v = vgALL, g = gstatALL,
                              model = model, fit.method=6, correct.diagonal = 1.01)
plot(vgALL, model=fitCrossed)

# Caso univariado: Podría hacerse de otras formas, pero para que el modelo pase en automático mediante
# las funciones proporcionadas el formato de salida debe ser una lista con nombres var.PC# y objetos
# de variogramModel (Por lo menos en lo que se refiere al uso de gstat)

model.PM10 <- list(
  'PM10.PC1' = vgm(psill = 20000, "Exp", range = 2103.25, nugget =  0),
  'PM10.PC2' = vgm(psill = 101494.96, "Exp", range = 1484.57, nugget = 0)
)


fitSingle <- lapply(names(gstatPM10$data), function(i) {
        gstat::fit.variogram(
          vgPM10[vgPM10$id == i,],
          model = model.PM10[[i]],
          fit.method = 6)
})

names(fitSingle) <- names(model.PM10)

plot(vgPM10[vgPM10$id == 'PM10.PC1', ], model = fitSingle$PM10.PC1)

# --------------------------------------------------------------------------------------------------
#                               PASO 2: Agregar los modelos a SpatFD
# --------------------------------------------------------------------------------------------------
# Básicamente agregar el modelo se da en dos pasos:
#  1. Llevar los modelos al lenguaje usual (Dataframe, parecido a gstat)
#  2. Generación de las funciones representantes de cada modelo (i.e. Combinación de funciones)
#
# Lo realmente importante es hacer que el dataframe esté bien para que la selección/combinación de
# modelos se solucione de forma correcta.
