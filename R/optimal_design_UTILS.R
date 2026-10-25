# ===================== FUNCIONES INTERNAS DE optimal_design() =================
# Solo las usa optimal_design() (y sus métodos). La explicación de lo que se
# minimiza y de las decisiones tomadas está en dev/03_muestreo_optimo.pdf.

# CANDIDATOS AL AZAR: Dentro del mapa de los puntos objetivo si lo hay; si no,
# dentro del rectángulo que contiene los puntos objetivo y las estaciones
random_candidates <- function(size, map, targets, stations){
  if (!is.null(map)) {
    area   <- sf::st_union(sf::st_geometry(map))
    points <- sf::st_coordinates(sf::st_sample(area, size = size, type = 'random'))[, 1:2, drop = F]
  } else {
    box    <- apply(rbind(targets, stations), 2, range)
    points <- sapply(seq_len(ncol(box)), function(j) stats::runif(size, box[1, j], box[2, j]))
  }
  dimnames(points) <- list(NULL, colnames(stations))
  points
}

# VARIANZA TOTAL DE UN DISEÑO: La suma, sobre los puntos objetivo, de la
# varianza de predicción de la curva con las estaciones existentes más las
# nuevas. Se usan los sistemas del modelo (los de 'scores' o los de 'lambda'):
# a cada campo se le agregan las estaciones nuevas y se arma de nuevo su
# sistema. La varianza no depende de los valores observados, solo de dónde
# están las estaciones, así que no hace falta inventar datos.
#
# Un diseño cuyo sistema no se puede resolver (por ejemplo, una estación nueva
# casi encima de otra sin efecto pepita) vale Inf: nunca se escoge.
design_variance <- function(systems, new, targets){
  total <- 0
  for (system in systems) {
    sites <- lapply(system$sites, function(s) rbind(s, new))
    built <- tryCatch(solver_system(system$fields, sites, system$variograms), error = function(e) NULL)
    if (is.null(built)) return(Inf)

    solution <- solver_solve(built, targets)
    total    <- total + sum(vapply(solution, function(field) sum(field$variance), numeric(1)))
  }
  total
}

# BARRA DE PROGRESO: Cuenta los candidatos probados de un total, con su
# etiqueta. Sin 'show' no muestra nada.
progress_bar <- function(label, total, show){
  if (!show) return(list(tick = function(n) invisible(NULL), close = function() invisible(NULL)))

  message(label)
  bar  <- utils::txtProgressBar(min = 0, max = total, style = 3)
  done <- 0
  list(tick  = function(n) {
         done <<- done + n
         utils::setTxtProgressBar(bar, done)
       },
       close = function() close(bar))
}

# ====================== VARIANZA RÁPIDA DE LOS DISEÑOS ========================
# design_variance() rearma y resuelve el sistema completo para cada diseño: es
# la referencia exacta, pero la búsqueda prueba muchísimos diseños. Aquí se usa
# que la matriz del sistema con estaciones nuevas es la del modelo con filas y
# columnas agregadas al final (bordeada):
#
#     A_S = | A   U |        b_S = | b |   (una columna por punto objetivo
#           | U'  D |              | c |    y campo que se predice)
#
#   A -> Sistema del modelo (ya factorizado), b -> sus lados derechos
#   U -> Variograma entre las observaciones existentes y las nuevas, más las
#        filas de las restricciones (las nuevas cuentan en su campo)
#   D -> Variograma entre las observaciones nuevas
#   c -> Variograma entre las observaciones nuevas y los puntos objetivo
#
# La varianza total es la suma de b' A^-1 b sobre los puntos objetivo. Con la
# inversa por bloques (complemento de Schur):
#
#     V(S) = V0 + traza( S_SS^-1 Q_SS ),   S = D - U' A^-1 U,   Q = R R',   R = c - U' A^-1 b
#
# donde V0 es la varianza sin estaciones nuevas. U, A^-1 U y R se calculan una
# sola vez para TODOS los candidatos. S y Q son de candidatos x candidatos, así
# que nunca se arman completas: solo las filas de las estaciones escogidas (y
# sus diagonales), que son las que hacen falta.
#
# AGREGAR UN CANDIDATO A LA VEZ: Para un diseño C y un candidato i, con la
# inversa por bloques de un elemento,
#
#     V(C + i) = V(C) + (Q_ii - 2 a' Q_Ci + a' Q_CC a) / s,   a = S_CC^-1 S_Ci,   s = S_ii - S_iC a
#
# y eso se calcula para TODOS los candidatos de una vez, con las filas de C. Es
# lo que usan la búsqueda voraz y los intercambios. Vale para grupos de un solo
# campo (lo común: cada score por separado, o la curva en 'lambda'); en grupos
# de varios campos (scores con variogramas cruzados) cada candidato se evalúa
# por separado, con matrices pequeñas.
#
# Devuelve dos funciones:
#   value(C) -> V(C) para los índices de los candidatos de un diseño
#   add(C)   -> V(C + i) para cada candidato i (Inf para los que ya están en C
#               o con los que el sistema no se puede resolver)
design_cache <- function(systems, candidates, targets, show = FALSE){
  N   <- nrow(candidates)
  bar <- progress_bar(paste('Preparing', N, 'candidates'), N * length(systems), show)
  on.exit(bar$close())

  parts <- lapply(systems, function(system) {
    fields <- system$fields
    m      <- length(fields)

    # Columnas del sistema para observar el campo j en 'points': variograma con
    # las observaciones existentes de cada campo y la fila de restricción de j
    columns <- function(points, j) {
      rbind(do.call(rbind, lapply(fields, function(i) {
              variogram_block(system$variograms[[i, j]], system$sites[[i]], points)
            })),
            matrix(as.numeric(fields == j), m, nrow(points)))
    }
    # Variogramas entre dos conjuntos de puntos, para cada pareja de campos
    pairs <- function(a, b) {
      do.call(rbind, lapply(fields, function(i) {
        do.call(cbind, lapply(fields, function(j) variogram_block(system$variograms[[i, j]], a, b)))
      }))
    }

    # Observaciones (filas de S y Q) de un conjunto de candidatos, en cada campo
    obs <- function(C) as.vector(outer(C, (seq_len(m) - 1) * N, '+'))

    U <- do.call(cbind, lapply(fields, columns, points = candidates))      # (n + m) x (N m)
    B <- do.call(cbind, lapply(fields, columns, points = targets))         # (n + m) x (T m)
    X <- qr.coef(system$qr, B)                                             # A^-1 b
    Y <- qr.coef(system$qr, U)                                             # A^-1 U

    # R, (N m) x (T m), por bloques de 1000 candidatos: así la memoria no se
    # multiplica con las copias temporales de los variogramas
    R <- matrix(0, N * m, ncol(B))
    for (rows in split(seq_len(N), ceiling(seq_len(N) / 1000))) {
      R[obs(rows), ] <- pairs(candidates[rows, , drop = F], targets) - crossprod(U[, obs(rows), drop = F], X)
      bar$tick(length(rows))
    }

    V0 <- sum(B * X)
    dS <- -colSums(U * Y)                       # Diagonal de S (el variograma vale 0 en el origen)
    dQ <- rowSums(R^2)                          # Diagonal de Q

    # Filas de S y Q de los candidatos C
    S_rows <- function(C) pairs(candidates[C, , drop = F], candidates) - crossprod(U[, obs(C), drop = F], Y)
    Q_rows <- function(C) tcrossprod(R[obs(C), , drop = F], R)

    # V(C) en este grupo
    value <- function(C) {
      if (length(C) == 0) return(V0)
      S   <- S_rows(C)[, obs(C), drop = F]
      inv <- tryCatch(solve(S), error = function(e) NULL)
      if (is.null(inv)) return(Inf)
      V0 + sum(inv * Q_rows(C)[, obs(C), drop = F])
    }

    # V(C + i) en este grupo, para todos los candidatos i
    add <- function(C) {
      if (m > 1) {
        return(vapply(seq_len(N), function(i) if (i %in% C) Inf else value(c(C, i)), numeric(1)))
      }

      if (length(C) == 0) {
        s   <- dS
        out <- V0 + dQ / s
      } else {
        SC  <- S_rows(C)
        QC  <- Q_rows(C)
        inv <- tryCatch(solve(SC[, C, drop = F]), error = function(e) NULL)
        if (is.null(inv)) return(rep(Inf, N))

        a   <- inv %*% SC                                         # Una columna a por candidato
        s   <- dS - colSums(SC * a)
        num <- dQ - 2 * colSums(a * QC) + colSums(a * (QC[, C, drop = F] %*% a))
        out <- V0 + sum(inv * QC[, C, drop = F]) + num / s
      }

      # Agregar una estación solo puede bajar la varianza: s < 0. Un s casi cero
      # es un candidato que no se puede agregar (por ejemplo, encima de otro)
      out[s > -1e-10 * max(abs(dS))] <- Inf
      out[C] <- Inf
      out
    }

    list(value = value, add = add)
  })

  list(value = function(C) sum(vapply(parts, function(part) part$value(C), numeric(1))),
       add   = function(C) Reduce(`+`, lapply(parts, function(part) part$add(C))))
}

# BÚSQUEDA VORAZ: Agrega las estaciones una por una; cada vez escoge, entre los
# candidatos que quedan, el que más reduce la varianza total. Devuelve los
# índices de los candidatos escogidos y la varianza después de cada uno.
greedy_design <- function(cache, k, N, show = FALSE){
  chosen   <- integer(0)
  variance <- numeric(0)
  bar      <- progress_bar(paste('Adding', k, 'stations: candidates tried'), k * N, show)
  on.exit(bar$close())

  for (step in seq_len(k)) {
    values   <- cache$add(chosen)
    chosen   <- c(chosen, which.min(values))
    variance <- c(variance, min(values))
    bar$tick(N)
  }

  list(chosen = chosen, variance = variance)
}

# INTERCAMBIOS: Para cada estación nueva, busca el candidato que más baja la
# varianza si la reemplaza, y la cambia si mejora. Repite mientras haya mejoras
# (o hasta 'sweeps' vueltas). La varianza nunca sube.
exchange_design <- function(cache, chosen, best, N, sweeps = 10, show = FALSE){
  for (sweep in seq_len(sweeps)) {
    improved <- FALSE
    bar      <- progress_bar(paste0('Exchanges, round ', sweep, ': candidates tried'), length(chosen) * N, show)

    for (j in seq_along(chosen)) {
      values <- cache$add(chosen[-j])                 # Reemplazar la estación j por cada candidato
      values[chosen] <- Inf
      if (min(values) < best * (1 - 1e-10)) {
        chosen[j] <- which.min(values)
        best      <- min(values)
        improved  <- TRUE
      }
      bar$tick(N)
    }

    bar$close()
    if (!improved) break
  }

  list(chosen = chosen, variance = best)
}
