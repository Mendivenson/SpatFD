# ============================ SOLVER DE KRIGING ===============================
# Funciones internas que resuelven los sistemas de kriging y cokriging, simple u
# ordinario, a partir de los modelos compilados. Las usan predict() y
# leave_one_out().
#
# KRIGING SIMPLE: la media es conocida e igual a cero (los scores están
# centrados). KRIGING ORDINARIO: la media es una constante desconocida; se
# agrega una restricción por campo (o por variable, en lambda): los pesos del
# campo que se predice suman 1 y los de cada uno de los demás suman 0.
#
# Con esas restricciones, sumar una constante a un bloque de covarianzas no
# cambia ni los pesos ni la varianza: la constante se cancela con los
# multiplicadores de Lagrange. Por eso el kriging ordinario puede usar modelos
# sin meseta, con la pseudo-covarianza C(h) = -gamma(h).

# Modelos sin meseta: solo se pueden usar con kriging ordinario
unbounded_models <- c('linear', 'power', 'logarithmic', 'spline', 'legendre')

# COVARIANZAS DE LOS SCORES: A partir de un modelo ya compilado, ver
# add_model(), se construye una matriz-lista con la función de covarianza
# de cada pareja de campos, C(h) = meseta - gamma(h). Las parejas sin
# variograma cruzado quedan en NULL, es decir, se toman como independientes.
# Con kriging ordinario, un variograma sin meseta usa la pseudo-covarianza
# C(h) = -gamma(h).
score_covariances <- function(compiled, ids, ordinary = FALSE){
  q    <- length(ids)
  covs <- matrix(vector('list', q * q), q, q, dimnames = list(ids, ids))

  for (i in seq_len(q)) {
    for (j in i:q) {
      # El id de un variograma cruzado puede venir en cualquiera de los dos órdenes
      key <- if (i == j) ids[i] else paste(ids[c(i, j)], ids[c(j, i)], sep = '.')
      key <- intersect(key, names(compiled$variograms))
      if (length(key) == 0) next

      structures <- compiled$table$model[compiled$table$var == key[1]]
      bounded    <- !any(structures %in% unbounded_models)
      if (!bounded && !ordinary) {
        stop("The model of '", key[1], "' has no sill (",
             paste(intersect(structures, unbounded_models), collapse = ', '),
             "), so it cannot be used for simple kriging. Use kriging = 'ordinary'", call. = FALSE)
      }
      if (bounded && is.na(compiled$sills[[key[1]]])) {
        stop("The model of '", key[1], "' has missing values in 'sigma'", call. = FALSE)
      }

      covs[[i, j]] <- covs[[j, i]] <- local({
        sill  <- if (bounded) compiled$sills[[key[1]]] else 0
        gamma <- compiled$variograms[[key[1]]]
        function(h) sill - gamma(h)
      })
    }
  }
  covs
}

# DISTANCIAS ENTRE DOS CONJUNTOS DE SITIOS: Se calcula coordenada a coordenada
# para que la distancia entre dos sitios iguales sea exactamente cero (el
# efecto pepita depende de ello).
cross_dist <- function(a, b){
  d2 <- matrix(0, nrow(a), nrow(b))
  for (k in seq_len(ncol(a))) {
    d2 <- d2 + outer(a[, k], b[, k], '-')^2
  }
  sqrt(d2)
}

# SISTEMA DE KRIGING CON RESTRICCIONES: Resuelve
#
#     | Sigma  F |  | W |   | C0 |
#     | F'     0 |  | M | = | E  |
#
# donde F tiene una columna por campo (1 en las filas de ese campo) y E indica
# el campo objetivo de cada columna de C0. Devuelve los pesos W y la varianza
# de predicción, C(0) - suma por columnas de [C0; E] * [W; M]. La matriz ya no
# es definida positiva, así que se resuelve por LU.
constrained_solve <- function(Sigma, C0, rows, cols, sills){
  k <- length(rows)
  F <- matrix(0, nrow(Sigma), k)
  E <- matrix(0, k, ncol(C0))
  for (a in seq_len(k)) {
    F[rows[[a]], a] <- 1
    E[a, cols[[a]]] <- 1
  }

  system <- rbind(cbind(Sigma, F), cbind(t(F), matrix(0, k, k)))
  rhs    <- rbind(C0, E)
  sol    <- tryCatch(solve(system, rhs), error = function(e) {
    stop('The ordinary kriging system is singular. ',
         'Check for duplicated sites or for an invalid model', call. = FALSE)
  })

  list(weights  = sol[seq_len(nrow(Sigma)), , drop = FALSE],
       variance = sills - colSums(rhs * sol))
}

# SOLVER DE KRIGING DE CAMPOS ESCALARES (los scores):
#   sites     -> Lista con las coordenadas de los sitios de cada campo
#   values    -> Lista con los valores observados de cada campo
#   covs      -> Matriz-lista de funciones de covarianza, ver score_covariances()
#   newcoords -> Coordenadas de los sitios de predicción
#   ordinary  -> Kriging ordinario (media constante desconocida) o simple
#                (media cero)
#
# Sigma se arma por bloques, uno por pareja de campos. Sin covarianzas cruzadas
# es diagonal por bloques y el sistema equivale a un kriging por campo; con
# ellas es el sistema de cokriging. Los casos univariado y multivariado son
# entonces el mismo sistema.
#
# Los campos que no están conectados por ninguna covarianza cruzada no
# comparten información, así que cada grupo de campos conectados se resuelve
# por aparte: da lo mismo que resolver el sistema completo y evita factorizar
# una matriz llena de ceros.
solve_kriging <- function(sites, values, covs, newcoords, ordinary = FALSE){
  q <- length(sites)
  m <- nrow(newcoords)

  # Los modelos no siempre conservan las dimensiones de la matriz de distancias.
  # Si las distancias son entre un conjunto de sitios y él mismo la matriz es
  # simétrica y basta evaluar la mitad.
  evaluate <- function(covariance, d, symmetric = FALSE){
    if (!symmetric) return(matrix(covariance(d), nrow(d), ncol(d)))

    half <- lower.tri(d, diag = TRUE)
    out  <- matrix(0, nrow(d), ncol(d))
    out[half] <- covariance(d[half])
    out[upper.tri(out)] <- t(out)[upper.tri(out)]
    out
  }

  # DISTANCIAS: Los scores de una misma variable comparten sitios, así que cada
  # matriz de distancias se calcula una sola vez por conjunto distinto de sitios
  places <- list()
  place  <- integer(q)                   # Conjunto de sitios de cada campo
  for (i in seq_len(q)) {
    known <- Position(function(p) identical(p, sites[[i]]), places)
    if (is.na(known)) {
      places <- c(places, list(sites[[i]]))
      known  <- length(places)
    }
    place[i] <- known
  }

  computed <- list()
  distance <- function(a, b){            # b = 0 son los sitios de predicción
    key <- paste(a, b)
    if (is.null(computed[[key]])) {
      computed[[key]] <<- cross_dist(places[[a]], if (b == 0) newcoords else places[[b]])
    }
    computed[[key]]
  }

  # GRUPOS: Campos conectados, directa o indirectamente, por covarianzas cruzadas
  group <- seq_len(q)
  for (i in seq_len(q)) {
    for (j in seq_len(q)) {
      if (!is.null(covs[[i, j]])) group[group == group[j]] <- group[i]
    }
  }

  pred <- variance <- matrix(0, m, q)

  for (g in unique(group)) {
    fields <- which(group == g)
    k <- length(fields)
    n <- vapply(sites[fields], nrow, integer(1))

    rows <- split(seq_len(sum(n)), rep(seq_len(k), n))         # Filas de cada campo
    cols <- split(seq_len(m * k), rep(seq_len(k), each = m))   # Columnas de cada campo objetivo

    Sigma <- matrix(0, sum(n), sum(n))    # Covarianza entre observaciones
    C0    <- matrix(0, sum(n), m * k)     # Covarianza entre observaciones y sitios nuevos

    for (a in seq_len(k)) {
      for (b in a:k) {
        i <- fields[a]
        j <- fields[b]
        if (is.null(covs[[i, j]])) next
        same_sites <- place[i] == place[j]

        # Sigma es simétrica: solo se evalúa un bloque de cada pareja de campos
        block <- evaluate(covs[[i, j]], distance(place[i], place[j]), symmetric = same_sites)
        Sigma[rows[[a]], rows[[b]]] <- block
        if (a < b) Sigma[rows[[b]], rows[[a]]] <- t(block)

        # La covarianza con los sitios nuevos solo depende de los sitios del
        # campo observado: si ambos campos comparten sitios, es la misma
        block <- evaluate(covs[[i, j]], distance(place[i], 0))
        C0[rows[[a]], cols[[b]]] <- block
        if (a < b) {
          C0[rows[[b]], cols[[a]]] <- if (same_sites) block else evaluate(covs[[i, j]], distance(place[j], 0))
        }
      }
    }

    sills <- vapply(fields, function(j) covs[[j, j]](0), numeric(1))

    # KRIGING ORDINARIO: sistema con una restricción por campo del grupo
    if (ordinary) {
      solved <- constrained_solve(Sigma, C0, rows, cols, rep(sills, each = m))
      pred[, fields]     <- matrix(crossprod(solved$weights, unlist(values[fields])), m, k)
      variance[, fields] <- matrix(solved$variance, m, k)
      next
    }

    # KRIGING SIMPLE: Una sola factorización sirve para todos los sitios y
    # campos del grupo. Cholesky falla si Sigma no es definida positiva, lo que
    # delata un modelo inválido.
    R <- tryCatch(chol(Sigma), error = function(e) {
      stop('The covariance matrix of the scores is not positive definite. ',
           'Check for duplicated sites or for an invalid model', call. = FALSE)
    })
    #
    # Con Sigma = R'R y V = R'^-1 C0 no hace falta calcular los pesos Sigma^-1 C0:
    #   predicción = C0' Sigma^-1 z  = V' (R'^-1 z)
    #   varianza   = C(0) - C0' Sigma^-1 C0 = C(0) - suma por columnas de V^2
    V <- backsolve(R, C0, transpose = TRUE)
    z <- backsolve(R, unlist(values[fields]), transpose = TRUE)

    pred[, fields]     <- matrix(crossprod(V, z), m, k)
    variance[, fields] <- matrix(rep(sills, each = m) - colSums(V^2), m, k)
  }

  list(pred = pred, variance = variance)
}
