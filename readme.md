# SpatFD v2

[![CRAN status](https://www.r-pkg.org/badges/version/SpatFD)](https://cran.r-project.org/package=SpatFD)
[![Licence](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Estado](https://img.shields.io/badge/estado-desarrollo%20activo%20(fase%20temprana)-orange.svg)]()

> ⚠️ **En desarrollo activo, fase temprana.** Este repositorio es un fork/rediseño en desarrollo de la v2 de SpatFD y todavía no debe usarse en producción.
>
> - Versión estable actual (CRAN): [SpatFD en CRAN](https://cran.r-project.org/web/packages/SpatFD/index.html)
> - Repositorio original del paquete: [mpbohorquezc/SpatFD-Functional-Geostatistics](https://github.com/mpbohorquezc/SpatFD-Functional-Geostatistics)

## Descripción

Este repositorio contiene el desarrollo de la segunda versión de **SpatFD**, paquete de R para geoestadística funcional (predicción espacial de curvas mediante *functional kriging* y *cokriging*, muestreo óptimo, clasificación supervisada y simulación de datos funcionales espacialmente correlacionados).

El objetivo principal de esta v2 es **estandarizar el desarrollo** a lo largo de todas las utilidades del paquete, haciéndolo:

- **Más mantenible**, con una arquitectura interna consistente.
- **Más amigable** para personas sin conocimiento extenso de geoestadística.

Los conceptos y algoritmos de modelamiento geoestadístico funcional no recibirán cambios importantes. El cambio principal es **eliminar a `gstat` como columna vertebral del paquete**, aunque se mantendrá la posibilidad de usarlo como una opción más.

## Funcionalidades

Se conservarán las funcionalidades del paquete original, pero se separa la responsabilidad del modelamiento:

- El paquete se encargará únicamente del **manejo de datos funcionales** (construcción de bases funcionales, suavizado, análisis funcional, etc.).
- El **modelamiento geoestadístico** se desacopla del paquete para dar mayor flexibilidad al usuario.
- La idea es que el usuario pueda conectar el output funcional con `gstat`, `geoR`, o incluso con **modelos personalizados**.

## Roadmap / TODO

### 23/Agosto

- ~~Hacer el espejo de los modelos de gstat en modelos de SpatFD.~~
- ~~Generar una función para listar y graficar todos los modelos.~~
- ~~Generar las funciones de gstat a y la del switch para el caso de clases simples como listas~~
- ~~Dejar un test pipeline listo para predicción.~~ 

### 30/Agosto

- Modificación comportamiento operador `<-`
- Predicción univariada
- Predicción multivariada

### 6/Septiembre

- Prediction methods output
- Documentación 
   - data inclussion 
   - data examples inclussion
- vignettes docs 

### 13/Septiembre

- Preparación exposición SpatFD - Encuentro de semilleros


### 20/Septiembre 

- Cross validation test

### 27/Septiembre 

- Muestreo óptimo 

> Este roadmap corresponde a una fase temprana del desarrollo y cambiará conforme avance el proyecto.

## Desarrollo

Este proyecto se está desarrollando en conjunto con la autora original del paquete ([mpbohorquezc](https://github.com/mpbohorquezc)). Una vez finalizado, esta v2 reemplazará la versión actualmente publicada en CRAN.

## Referencias

El marco teórico de geoestadística funcional en el que se basa este paquete:

- Bohórquez, M., Giraldo, R., & Mateu, J. (2016). Optimal sampling for spatial prediction of functional data. *Statistical Methods & Applications*, 25(1), 39-54.
- Bohórquez, M., Giraldo, R., & Mateu, J. (2017). Multivariate functional random fields: prediction and optimal sampling. *Stochastic Environmental Research and Risk Assessment*, 31, 53–70.
- Bohórquez M., Giraldo R. and Mateu J. (2021). *Geostatistical Functional Data Analysis: Theory and Methods*. John Wiley & Sons, Chichester, UK. ISBN: 978-1-119-38784-8.

## Licencia

Este paquete es software libre y de código abierto, licenciado bajo [GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html).