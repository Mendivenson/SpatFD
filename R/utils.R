# ======================= FUNCIONES INTERNAS COMPARTIDAS =======================
# Todo lo funcional del objeto (datos originales, curvas suavizadas y
# componentes principales) vive en un único atributo oculto, 'internal', que es
# una lista normal:
#
#   coords -> Coordenadas de los sitios con datos
#   data   -> Matriz de datos (una fila por tiempo, una columna por sitio)
#   time   -> Tiempos tal como los dio el usuario (números o fechas)
#   t      -> Tiempos numéricos usados en el suavizamiento, ver spatfd_UTILS.R
#   fd     -> Curvas suavizadas, salida de fda
#   pca    -> Componentes principales, salida de fda::pca.fd()
#   args   -> Argumentos del suavizamiento, para poder repetirlo
#
# El resto del paquete accede a él solo a través de estas funciones, así que
# el usuario nunca tiene que manipular objetos del paquete fda.

# CONSTRUCTOR: Arma el objeto a partir de sus partes. El nombre de la variable
# queda en su propio atributo, visible, para nombrar gráficos y modelos.
new_spatfd <- function(scores, variance, missing, internal, name){
  structure(list(scores = scores, variance = variance, missing = missing),
            class    = 'SpatFD',
            name     = name,
            internal = internal)
}

# ACCESORES
get_name     <- function(sfd) attr(sfd, 'name')
get_internal <- function(sfd, what) attr(sfd, 'internal')[[what]]
get_coords   <- function(sfd) get_internal(sfd, 'coords')
get_fd       <- function(sfd) get_internal(sfd, 'fd')
get_pca      <- function(sfd) get_internal(sfd, 'pca')

# TIEMPOS EN DÍAS: Las fechas se pasan a días desde 1970-01-01, con la hora
# como fracción del día (como en Excel, con otro origen). Siempre en UTC: así
# los tiempos crecen siempre, sin horas repetidas ni saltadas por el cambio de
# horario, y un objeto con fechas y otro con fecha y hora quedan en el mismo
# eje (una fecha es las 00:00 UTC de ese día). Es solo para los cálculos: el
# usuario sigue viendo sus tiempos tal cual. Los números se usan tal cual.
as_days <- function(time){
  if (inherits(time, 'Date'))   return(as.numeric(time))
  if (inherits(time, 'POSIXt')) return(as.numeric(as.POSIXct(time)) / 86400)
  time
}

# DE DÍAS A TIEMPOS DEL USUARIO: La inversa de as_days(), con la clase (y la
# zona horaria) de los tiempos de 'template'
from_days <- function(days, template){
  if (inherits(template, 'Date'))   return(as.Date(days, origin = '1970-01-01'))
  if (inherits(template, 'POSIXt')) {
    tz <- attr(as.POSIXct(template), 'tzone')
    return(as.POSIXct(days * 86400, origin = '1970-01-01', tz = if (is.null(tz)) '' else tz[1]))
  }
  days
}

# CELDAS DE UNA GRILLA REGULAR: Los sitios de una predicción sobre un mapa
# forman una grilla regular. Para dibujarla como imagen se necesitan sus
# coordenadas x e y distintas y una matriz con un valor por celda (filas = y,
# columnas = x), con NA fuera del mapa. 'fill' pasa los valores de los sitios a
# esa matriz.
grid_cells <- function(grid){
  step <- function(v) {
    values <- sort(unique(signif(v, 8)))
    if (length(values) > 1) min(diff(values)) else 1
  }
  dx <- step(grid[, 1])
  dy <- step(grid[, 2])
  ix <- round((grid[, 1] - min(grid[, 1])) / dx) + 1
  iy <- round((grid[, 2] - min(grid[, 2])) / dy) + 1

  list(x    = min(grid[, 1]) + dx * (seq_len(max(ix)) - 1),
       y    = min(grid[, 2]) + dy * (seq_len(max(iy)) - 1),
       fill = function(values) {
         z <- matrix(NA_real_, max(iy), max(ix))
         z[cbind(iy, ix)] <- values
         z
       })
}

# VARIOS OBJETOS SpatFD: Revisa los objetos dados en '...' (por to_gstat() y
# spatfd_model()) y los devuelve en una lista nombrada con el nombre de cada
# variable. Los nombres deben ser distintos: forman los ids de los scores.
spatfd_list <- function(...){
  objects <- list(...)

  stopifnot(
    "At least one SpatFD object must be given" = length(objects) > 0,
    "All the objects in '...' must be of class SpatFD" = all(vapply(objects, inherits, logical(1), what = 'SpatFD'))
  )

  names(objects) <- vapply(objects, get_name, character(1))
  repeated <- unique(names(objects)[duplicated(names(objects))])
  if (length(repeated) > 0) {
    stop("Several objects are named ", paste0("'", repeated, "'", collapse = ', '),
         ". Give them different names with the argument 'name' of spatfd()", call. = F)
  }

  # Todas las variables deben estar en el mismo sistema de coordenadas
  coord_names <- lapply(objects, function(sfd) colnames(get_coords(sfd)))
  if (!all(vapply(coord_names, identical, logical(1), coord_names[[1]]))) {
    stop("The coordinates of all the objects must have the same column names", call. = F)
  }

  objects
}

# PRODUCTOS INTERNOS ENTRE CURVAS: Matriz con <f1_i, f2_j> integrado sobre
# 'rng'. Si las dos curvas están en la misma base y se integra en todo su
# intervalo, el resultado es exacto con la matriz de productos internos de la
# base. Si no, se integra por cuadratura de Gauss-Legendre (5 puntos en 'pieces'
# tramos), mucho más precisa que fda::inprod(), que se detiene en un error
# relativo de 1e-4.
inner_product <- function(fd1, fd2, rng, pieces = 1000){
  if (identical(fd1$basis, fd2$basis) && all(fd1$basis$rangeval == rng)) {
    return(crossprod(fd1$coefs, fda::eval.penalty(fd1$basis, fda::int2Lfd(0)) %*% fd2$coefs))
  }

  nodes   <- c(-0.9061798459386640, -0.5384693101056831, 0, 0.5384693101056831, 0.9061798459386640)
  weights <- c(0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891)

  breaks <- seq(rng[1], rng[2], length.out = pieces + 1)
  half   <- diff(breaks) / 2
  t      <- as.vector(outer(nodes, half) + rep(breaks[-1] - half, each = length(nodes)))
  w      <- as.vector(outer(weights, half))

  crossprod(fda::eval.fd(t, fd1) * w, fda::eval.fd(t, fd2))
}
