# ================================== DATOS =====================================
# Conjuntos de datos heredados de la versión CRAN del paquete. Cada archivo de
# data/ contiene varios objetos y se carga con data(): no se usa LazyData.

#' Air quality in Mexico City: PM10 and NO2
#'
#' Concentrations of particulate matter (PM10) and nitrogen dioxide (NO2)
#' measured over time at monitoring stations of Mexico City, with the
#' coordinates of the stations. The two pollutants are not measured at the
#' same set of stations, which makes the data suitable for functional
#' cokriging.
#'
#' @format `data(air_mexico)` loads five objects:
#' \describe{
#'   \item{`pm10_mexico`}{Numeric matrix with 4344 rows and 13 columns: PM10
#'     measurements, one row per measurement time and one column per station.}
#'   \item{`coords_pm10_mexico`}{Data frame with 13 rows and the columns `X` and `Y`:
#'     coordinates of the PM10 stations. Row names are the station codes, the
#'     column names of `pm10_mexico`.}
#'   \item{`no2_mexico`}{Numeric matrix with 4292 rows and 18 columns: NO2
#'     measurements, one row per measurement time and one column per station.}
#'   \item{`coords_no2_mexico`}{Data frame with 18 rows and the columns `X` and `Y`:
#'     coordinates of the NO2 stations. Row names are the station codes, the
#'     column names of `no2_mexico`.}
#'   \item{`map_mexico`}{Map of the municipalities of the area, an `sf` object
#'     with 54 multipolygons. The \pkg{sf} package is needed to work with it.}
#' }
#'
#' @examples
#' data(air_mexico)
#'
#' pm10 <- spatfd(pm10_mexico, coords = as.matrix(coords_pm10_mexico),
#'                basis = "fourier", nbasis = 21, lambda = 0.000001,
#'                nharm = 2, name = "PM10")
#' no2  <- spatfd(no2_mexico, coords = as.matrix(coords_no2_mexico),
#'                basis = "fourier", nbasis = 21, lambda = 0.000001,
#'                nharm = 2, name = "NO2")
#'
#' summary(pm10 + no2)
#'
#' @docType data
#' @keywords datasets
#' @name air_mexico
#' @aliases pm10_mexico coords_pm10_mexico no2_mexico coords_no2_mexico map_mexico
#' @usage data(air_mexico)
NULL

#' Air quality in Bogota, Colombia: PM10
#'
#' Particulate matter (PM10) measured over time at 10 monitoring stations of
#' Bogota, with the coordinates of the stations.
#'
#' @format `data(air_bogota)` loads three objects:
#' \describe{
#'   \item{`pm10_bogota`}{Data frame with 8761 rows and 10 columns: PM10
#'     measurements, one row per measurement time and one column per station.}
#'   \item{`coords_bogota`}{Data frame with 10 rows and the columns `ESTACION`
#'     (station), `X` and `Y`: coordinates of the stations, in the same order
#'     as the columns of `pm10_bogota`. The station names in `ESTACION` are
#'     abbreviated and do not match the column names of `pm10_bogota`.}
#'   \item{`map_bogota`}{Map of Bogota, an `sf` object with one polygon. The \pkg{sf}
#'     package is needed to work with it.}
#' }
#'
#' @source [Secretaría de Ambiente de Bogotá](https://www.ambientebogota.gov.co/es/inicio)
#'
#' @section Notes:
#' For real-time information visit the
#' [air quality monitoring network of Bogota](http://rmcab.ambientebogota.gov.co/home/map).
#'
#' @examples
#' data(air_bogota)
#'
#' # spatfd() matches the columns of the data with the row names of the coordinates
#' coords <- as.matrix(coords_bogota[, c("X", "Y")])
#' rownames(coords) <- colnames(pm10_bogota)
#'
#' sfd <- spatfd(as.matrix(pm10_bogota), coords = coords,
#'               basis = "bsplines", nbasis = 91, nharm = 3, name = "PM10")
#' summary(sfd)
#'
#' @docType data
#' @keywords datasets
#' @name air_bogota
#' @aliases pm10_bogota coords_bogota map_bogota
#' @usage data(air_bogota)
NULL

#' Vowels: spatially correlated functional data with classes
#'
#' Signals recorded at 21 sites arranged on a regular grid, labelled with one
#' of five classes. In the CRAN version of the package they are used for the
#' classification of spatially correlated functional data, with the classes
#' named after the vowels `a`, `e`, `i`, `o` and `u` and 228 repetitions per
#' class. The current version has no classification functions yet.
#'
#' @format `data(vowels)` loads two objects:
#' \describe{
#'   \item{`vowels`}{Numeric matrix with 17100 rows and 22 columns. The
#'     columns `E1` to `E21` are the 21 sites, and the last column, `CLASE`,
#'     is the class of each row, from 1 to 5, with 3420 rows per class.}
#'   \item{`coords_vowels`}{Data frame with 21 rows and the columns `x` and
#'     `y`: coordinates of the 21 sites, in the same order as the columns of
#'     `vowels`.}
#' }
#'
#' @examples
#' data(vowels)
#'
#' dim(vowels)
#' table(vowels[, "CLASE"])
#' plot(coords_vowels, pch = 19, main = "Sites of the vowels data")
#'
#' @docType data
#' @keywords datasets
#' @name vowels
#' @aliases coords_vowels
#' @usage data(vowels)
NULL
