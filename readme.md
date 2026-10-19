# SpatFD v2

[![CRAN status](https://www.r-pkg.org/badges/version/SpatFD)](https://cran.r-project.org/package=SpatFD)
[![Licence](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Estado](https://img.shields.io/badge/estado-desarrollo%20activo-orange.svg)]()

> ⚠️ **En desarrollo activo.** Este repositorio es un fork/rediseño en desarrollo de la v2 de SpatFD (versión `0.0.1.9000`). La predicción espacial ya funciona de principio a fin, pero la interfaz puede cambiar y todavía no debe usarse en producción.
>
> - Versión estable actual (CRAN): [SpatFD en CRAN](https://cran.r-project.org/web/packages/SpatFD/index.html)
> - Repositorio original del paquete: [mpbohorquezc/SpatFD-Functional-Geostatistics](https://github.com/mpbohorquezc/SpatFD-Functional-Geostatistics)

## Descripción

Este repositorio contiene el desarrollo de la segunda versión de **SpatFD**, paquete de R para geoestadística funcional: predicción espacial de curvas mediante *functional kriging* y *cokriging*.

El objetivo principal de esta v2 es **estandarizar el desarrollo** a lo largo de todas las utilidades del paquete, haciéndolo:

- **Más mantenible**, con una arquitectura interna consistente.
- **Más amigable** para personas sin conocimiento extenso de geoestadística.

Los conceptos y algoritmos de modelamiento geoestadístico funcional no reciben cambios importantes. El cambio principal es **eliminar a `gstat` como columna vertebral del paquete**:

- El paquete se encarga del **manejo de los datos funcionales** (suavizamiento y componentes principales funcionales) y de la **predicción**, que se resuelve de forma matricial dentro del paquete.
- El **ajuste de los variogramas** de los *scores* se hace por fuera, con la herramienta que el usuario prefiera. Hoy el paquete exporta los *scores* a `gstat` y trae de vuelta los modelos ajustados con `gstat`.

## Estado actual

El flujo de predicción está completo y es un paquete instalable, con documentación de todas las funciones y 335 pruebas automáticas (`testthat`).

| Paso | Función | Qué hace |
|:-----|:--------|:---------|
| Objeto funcional | `spatfd()`, `+` | Suaviza las curvas y calcula sus componentes principales funcionales. Las variables se combinan con `+`; cada una conserva sus propios sitios. |
| Exploración | `summary()`, `plot()`, `var_names()` | Varianza explicada, curvas suavizadas y componentes principales. |
| Exportar los *scores* | `togstat()` | Deja los *scores* listos para estimar y ajustar variogramas con `gstat`. |
| Registrar modelos | `add_model()`, `toSpatFD()` | Traduce los modelos ajustados con `gstat`, los compila y los guarda dentro del objeto. Un objeto puede guardar varios modelos. |
| Modelos de variograma | `models()`, `availablemodels()` | 18 modelos espejo de los de `gstat`, que reproducen sus valores. |
| Predicción | `predict()` | Kriging (sin variogramas cruzados) o cokriging (con ellos), ordinario o simple, en sitios nuevos o sobre un mapa `sf`. Método `"scores"` (kriging de los *scores*) o `"lambda"` (combinación de las curvas observadas; con variogramas cruzados entre variables es cokriging funcional). |
| Resultados | `print()`, `plot()` | Curvas predichas, varianza de predicción y mapas de ambas. |
| Validación cruzada | `leave_one_out()` | Compara modelos y métodos dejando por fuera un sitio a la vez. |

Datos incluidos: `air_mexico` (PM10 y NO2 en Ciudad de México), `air_bogota` (PM10 en Bogotá) y `vowels`.

Por defecto la predicción es por kriging ordinario: la media es una constante desconocida, así que da lo mismo que las curvas estén centradas o no, y se pueden usar modelos sin meseta. Con `kriging = "simple"` la media es cero, como en la versión CRAN. Los resultados del método `"scores"` coinciden con el kriging simple y con el ordinario de `gstat`, y con la versión CRAN (`KS_scores_lambdas()`, `crossval_loo()`) en el caso simple, con diferencias relativas entre 10⁻¹⁶ y 10⁻¹³. El método `"lambda"` ya no sigue a la versión CRAN: usa las mesetas del modelo, y sus dos versiones coinciden con las fórmulas escritas con las covarianzas de `gstat`. La comparación está en [`dev/02_comparacion_gstat.R`](dev/02_comparacion_gstat.R).

### Limitaciones conocidas

- **Media constante.** El kriging ordinario estima una media constante; todavía no hay kriging universal (media con tendencia espacial).
- **Método `"lambda"` con varias variables**: solo combina las curvas en los tiempos que comparten, así que el `t` de `spatfd()` debe estar en una escala común. Con intervalos distintos predice sobre el común (en los datos de México, PM10 tiene 4344 instantes y NO2 4292, y se supone que empiezan a la misma hora).
- **Solo se traen modelos de `gstat`**, y la anisotropía se descarta en la traducción.
- **Vecindario global y distancias euclidianas**: sirve para coordenadas proyectadas y un número moderado de sitios.
- Todavía no se han migrado desde la versión CRAN el **muestreo óptimo**, la **clasificación** ni la **simulación**.

El detalle está en [`dev/brief_salvedades.pdf`](dev/brief_salvedades.pdf).

## Instalación

```r
# install.packages("remotes")
remotes::install_github("Mendivenson/SpatFD")
```

Para predecir sobre mapas se necesita además el paquete `sf`.

## Ejemplo de uso

Predicción de las curvas de PM10 y NO2 en Ciudad de México. Los dos contaminantes no se miden en las mismas estaciones (13 y 18).

```r
library(SpatFD)
library(gstat)
data(air_mexico)

# 1. Suavizamiento y componentes principales funcionales de cada variable
pm10 <- spatfd(pm10_mexico, coords = as.matrix(coords_pm10_mexico),
               basis = "fourier", nbasis = 21, lambda = 0.000001,
               nharm = 2, name = "PM10")
no2  <- spatfd(no2_mexico, coords = as.matrix(coords_no2_mexico),
               basis = "fourier", nbasis = 21, lambda = 0.000001,
               nharm = 2, name = "NO2")

aire <- pm10 + no2
summary(aire)
#> --------- PCA information ----------
#>        basis no. PC  variance
#> PM10 fourier      2 0.9564145
#> NO2  fourier      2 0.9658895
#> ------------------------------------

# 2. Variogramas de los scores, por fuera del paquete (aquí con gstat):
#    un modelo lineal de corregionalización para PM10 y NO2
scores <- togstat(aire)
vario  <- variogram(scores)

inicial <- vgm(647677.1, "Wav", 23317.05)
inicial <- vgm(127633,   "Wav",  9408.63, add.to = inicial)
inicial <- vgm(100,      "Nug",        0, add.to = inicial)
lmc <- fit.lmc(v = vario, g = scores, model = inicial,
               fit.method = 6, correct.diagonal = 1.01)

# 3. Registrar los modelos en el objeto (devuelve el objeto actualizado, como update())
aire <- aire |>
  add_model(lmc, name = "lmc") |>
  add_model(list(                                     # solo PM10, sin variogramas cruzados
    PM10.PC1 = vgm(200000, "Exp", 8000, nugget = 20000),
    PM10.PC2 = vgm(60000,  "Sph", 15000)), name = "pm10")

aire$models
#>      variables                     structures no. scores cross    sill
#> lmc  NO2, PM10                   wave, nugget          4  TRUE 1214940
#> pm10      PM10 nugget, exponential, spherical          2 FALSE  280000

# 4. Predicción en sitios nuevos: el modelo 'lmc' tiene variogramas cruzados -> cokriging.
#    Por defecto es kriging ordinario; con kriging = "simple", la media de los scores es cero.
sitios <- cbind(X = c(480000, 490000, 500000),
                Y = c(2140000, 2150000, 2160000))

pred <- predict(aire, sitios, model = "lmc")
pred
#> Predicted curves at 3 sites for: NO2, PM10
#> Method: scores, ordinary kriging, with cross-variograms
#>
#> Prediction variance:
#>           NO2     PM10
#> new1 136630.6 402843.9
#> new2 135366.3 395450.2
#> new3 169406.6 443848.3

pred$PM10$pred                  # curvas predichas: una fila por tiempo, una columna por sitio
plot(pred, variable = "PM10")   # predichas en color sobre las observadas en gris

# 5. Predicción sobre un mapa: plot() dibuja la predicción y su varianza lado a lado
mapa <- predict(aire, map_mexico, model = "lmc", variable = "PM10", n = 2000)
plot(mapa)
plot(mapa, time = 1000, main = "PM10")   # el mismo mapa en otro instante

# 6. Comparación de modelos y métodos por validación cruzada
leave_one_out(aire, variable = "PM10")
#> 13 sites, ordinary kriging, error = distance between the observed and the predicted curve
#>
#>   model method crossed     mean       sd      max
#> 1   lmc lambda    TRUE 486.1857 272.8846 1147.095
#> 2  pm10 scores   FALSE 556.2050 355.3138 1459.726
#> 3  pm10 lambda   FALSE 563.4934 353.2229 1467.241
#> 4   lmc scores    TRUE 595.5434 279.3465 1369.052
```

El modelo `pm10` usa parámetros fijos, solo como ilustración. Con el modelo `lmc`, el método `"lambda"` es cokriging funcional y solo cubre las horas comunes a PM10 y NO2, así que su error no es directamente comparable con el de los demás. Un recorrido más completo, paso a paso, está en [`dev/01_uso_paquete.R`](dev/01_uso_paquete.R).

## Roadmap / TODO

### Hecho

- ~~Espejo de los modelos de `gstat` en modelos de SpatFD, con función para listarlos y graficarlos.~~
- ~~Exportación de los *scores* a `gstat` y traducción de los modelos ajustados.~~
- ~~Estructura de paquete y registro de modelos dentro del objeto.~~
- ~~Predicción univariada (kriging) y multivariada (cokriging) con *solver* propio.~~
- ~~Kriging ordinario (por defecto), que admite modelos sin meseta.~~
- ~~Método `lambda` univariado, centrado, y cokriging funcional con varias variables.~~
- ~~Salida de la predicción: resumen, gráficos de curvas y mapas.~~
- ~~Datos y ejemplos incluidos en el paquete.~~
- ~~Validación cruzada dejando un sitio por fuera.~~

### Pendiente

- Traductores para modelos ajustados con otras herramientas (`geoR`, modelos personalizados).
- Muestreo óptimo.
- Clasificación y simulación de datos funcionales espacialmente correlacionados.
- Vignettes.

> Este roadmap cambiará conforme avance el proyecto.

## Desarrollo

Este proyecto se está desarrollando en conjunto con la autora original del paquete ([mpbohorquezc](https://github.com/mpbohorquezc)). Una vez finalizado, esta v2 reemplazará la versión actualmente publicada en CRAN.

```r
devtools::load_all()   # cargar el paquete desde el repositorio
devtools::test()       # correr las pruebas
```

La carpeta [`dev/`](dev) contiene los scripts de uso y de comparación con `gstat`, y los documentos técnicos de la implementación. Para entender cómo está construido el paquete (archivos, estructura del objeto, recorrido de una predicción y cómo extenderlo), el punto de partida es [`dev/estructura_paquete.pdf`](dev/estructura_paquete.pdf).

## Referencias

El marco teórico de geoestadística funcional en el que se basa este paquete:

- Bohórquez, M., Giraldo, R., & Mateu, J. (2016). Optimal sampling for spatial prediction of functional data. *Statistical Methods & Applications*, 25(1), 39-54.
- Bohórquez, M., Giraldo, R., & Mateu, J. (2017). Multivariate functional random fields: prediction and optimal sampling. *Stochastic Environmental Research and Risk Assessment*, 31, 53–70.
- Bohórquez M., Giraldo R. and Mateu J. (2021). *Geostatistical Functional Data Analysis: Theory and Methods*. John Wiley & Sons, Chichester, UK. ISBN: 978-1-119-38784-8.

## Licencia

Este paquete es software libre y de código abierto, licenciado bajo [GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html).
