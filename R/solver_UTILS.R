# ============================ SOLVER DE KRIGING ===============================
# Funciones internas (no son del usuario) que arman y resuelven el sistema de
# kriging ordinario. Las usan spatfd_model(), que arma, valida y factoriza el
# sistema una sola vez, y predict(), que lo resuelve para los sitios nuevos.
# La explicación completa del sistema está en dev/02_solver.pdf.
#
# Un CAMPO es lo que se predice: un score en el método 'scores' o la curva de
# una variable en el método 'lambda'. Cada campo tiene sus propios sitios
# observados. Para el solver da igual qué sea: solo necesita
#
#   sites      -> Lista con las coordenadas de los sitios de cada campo
#   variograms -> Matriz-lista (campos x campos) con el variograma de cada
#                 pareja de campos: una función de la distancia, o NULL si la
#                 pareja no tiene variograma (los campos son independientes)
#
# EL SISTEMA: Con n observaciones en total (los sitios de todos los campos,
# uno debajo del otro) y k campos,
#
#     | Gamma  F |   | W |   | G0 |
#     | F'     0 | x | M | = | E  |
#
#   Gamma (n x n) -> Variograma entre cada pareja de observaciones, por bloques
#   F     (n x k) -> Una columna por campo, con 1 en las filas de ese campo
#   G0            -> Variograma entre las observaciones y cada sitio nuevo
#   E             -> Qué campo se predice en cada columna de G0
#   W             -> Los pesos. M -> Los multiplicadores de Lagrange
#
# Las filas de arriba piden que el error sea mínimo; las de abajo, que los
# pesos del campo que se predice sumen 1 y los de los demás campos sumen 0
# (kriging ordinario). Con esas restricciones el sistema se puede escribir con
# los variogramas, sin la meseta.

# =============================== GRUPOS =======================================
# Campos que no están conectados por ningún variograma cruzado, ni directa ni
# indirectamente, no comparten información: cada grupo se resuelve aparte. Es
# lo mismo que resolver todo junto, pero sin matrices llenas de ceros.
solver_groups <- function(variograms){
  q     <- nrow(variograms)
  group <- seq_len(q)                           # Al inicio, cada campo es su propio grupo

  # Dos campos con variograma cruzado quedan en el mismo grupo
  for (i in seq_len(q)) {
    for (j in seq_len(q)) {
      if (!is.null(variograms[[i, j]])) group[group == group[j]] <- group[i]
    }
  }

  unname(split(rownames(variograms), group))
}

# ============================ DISTANCIAS ======================================
# Distancia euclidiana entre cada sitio de 'a' y cada sitio de 'b'
cross_dist <- function(a, b){
  d2 <- matrix(0, nrow(a), nrow(b))
  for (k in seq_len(ncol(a))) d2 <- d2 + outer(a[, k], b[, k], '-')^2
  sqrt(d2)
}

# BLOQUE DE VARIOGRAMAS: El variograma de la pareja de campos (i, j) entre los
# sitios 'a' y 'b'. Sin variograma, el bloque es cero.
variogram_block <- function(variogram, a, b){
  if (is.null(variogram)) return(matrix(0, nrow(a), nrow(b)))
  matrix(variogram(cross_dist(a, b)), nrow(a), nrow(b))
}

# ====================== ARMAR, VALIDAR Y FACTORIZAR ===========================
# Arma la matriz del sistema de un grupo de campos, revisa que el modelo sea
# válido y la factoriza (QR) para resolverla después las veces que haga falta.
# La matriz solo depende de los sitios observados y del modelo.
solver_system <- function(fields, sites, variograms){
  sites      <- sites[fields]
  variograms <- variograms[fields, fields, drop = F]

  n    <- vapply(sites, nrow, integer(1))
  rows <- split(seq_len(sum(n)), rep(seq_along(fields), n))    # Filas de cada campo
  names(rows) <- fields

  # GAMMA: Un bloque por pareja de campos
  Gamma <- matrix(0, sum(n), sum(n))
  for (i in fields) {
    for (j in fields) {
      Gamma[rows[[i]], rows[[j]]] <- variogram_block(variograms[[i, j]], sites[[i]], sites[[j]])
    }
  }

  # F: Las restricciones, una por campo
  F <- matrix(0, sum(n), length(fields))
  for (k in seq_along(fields)) F[rows[[k]], k] <- 1

  solver_check(Gamma, F)

  A <- rbind(cbind(Gamma, F),
             cbind(t(F), matrix(0, length(fields), length(fields))))

  list(fields = fields, sites = sites, variograms = variograms, rows = rows, qr = qr(A))
}

# VALIDEZ: El modelo es válido en los sitios observados si toda combinación de
# observaciones cuyos pesos suman cero en cada campo tiene varianza positiva:
# -Gamma debe ser definida positiva sobre esas combinaciones. Se revisan los
# valores propios de -Gamma proyectada sobre ellas (las columnas de Q, que son
# ortogonales a F).
#   Valor propio negativo  -> modelo no válido
#   Valor propio casi cero -> sistema casi singular
solver_check <- function(Gamma, F){
  if (nrow(Gamma) == ncol(F)) return(invisible(TRUE))          # Un sitio por campo: nada que revisar

  Q      <- qr.Q(qr(F), complete = TRUE)[, -seq_len(ncol(F)), drop = F]
  values <- eigen(-t(Q) %*% Gamma %*% Q, symmetric = TRUE, only.values = TRUE)$values
  scale  <- max(abs(values))

  if (min(values) < -1e-8 * scale) {
    stop('The model is not valid', call. = F)
  }
  if (min(values) <= 1e-12 * scale) {
    stop('The kriging system is singular or close to singular. Adding a nugget effect to the ',
         'model usually solves it', call. = F)
  }
  invisible(TRUE)
}

# ============================== RESOLVER ======================================
# Resuelve el sistema de un grupo para los sitios nuevos: una columna por cada
# pareja (campo que se predice, sitio nuevo). Devuelve, para cada campo que se
# predice, los pesos de todas las observaciones (filas) para cada sitio nuevo
# (columnas) y la varianza de predicción en cada sitio nuevo.
solver_solve <- function(system, newcoords){
  fields <- system$fields
  m      <- nrow(newcoords)
  n      <- length(unlist(system$rows))
  k      <- length(fields)

  out <- lapply(fields, function(target) {
    # LADO DERECHO: Variograma entre cada observación y cada sitio nuevo, más
    # las restricciones (1 para el campo que se predice, 0 para los demás)
    G0 <- do.call(rbind, lapply(fields, function(i) {
      variogram_block(system$variograms[[i, target]], system$sites[[i]], newcoords)
    }))
    E  <- matrix(as.numeric(fields == target), k, m)

    rhs      <- rbind(G0, E)
    solution <- qr.coef(system$qr, rhs)

    # VARIANZA: sum(pesos * variogramas) + multiplicador del campo que se predice
    list(weights  = solution[seq_len(n), , drop = F],
         variance = colSums(rhs * solution))
  })

  stats::setNames(out, fields)
}

# COMBINACIÓN: La predicción de un campo es la suma, campo por campo, de sus
# valores observados por sus pesos. 'values' tiene, para cada campo, una matriz
# con una columna por sitio observado: una fila para un score, una fila por
# tiempo para una curva.
solver_combine <- function(system, solution, values){
  Reduce(`+`, lapply(system$fields, function(i) {
    values[[i]] %*% solution$weights[system$rows[[i]], , drop = F]
  }))
}
