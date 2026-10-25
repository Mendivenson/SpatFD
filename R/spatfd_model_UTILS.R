# ==================== FUNCIONES INTERNAS DE spatfd_model() ====================
# Solo las usa spatfd_model() (y sus métodos). Lo compartido por varias
# funciones está en utils.R.
#
# Cada forma de dar el modelo (gstat, data frame o funciones) termina en la
# misma representación: una entrada por id de variograma con
#
#   variogram  -> Función de la distancia, conserva la forma de 'h' (vector o matriz)
#   sill       -> Meseta C(0); NA si el modelo no tiene meseta
#   structures -> Texto con las estructuras del modelo, para los resúmenes
#   source     -> 'gstat', 'data frame' o 'function'
#
# Así el solver no necesita saber de dónde vino cada modelo.

# ======================== NORMALIZACIÓN DEL MODELO ============================
#   model -> El modelo como lo dio el usuario
#   scale -> Tamaño del área de los sitios, para la meseta de las funciones
normalize_model <- function(model, scale){
  named <- function(x) !is.null(names(x)) && all(names(x) != '') && !anyDuplicated(names(x))

  if (inherits(model, 'gstat')) {
    if (length(model$model) == 0) stop('The gstat object has no fitted models', call. = F)
    return(from_gstat(model$model))
  }

  # Un variogramModel también es un data frame: se revisa antes
  if (inherits(model, 'variogramModel')) {
    stop('A single variogram model has no id. Give a named list of them, one per score ',
         '(for example list(PM10.PC1 = model))', call. = F)
  }
  if (is.data.frame(model)) return(from_table(model))

  if (is.list(model) && length(model) > 0 && named(model)) {
    if (all(vapply(model, inherits, logical(1), what = 'variogramModel'))) return(from_gstat(model))
    if (all(vapply(model, is.function, logical(1))))                       return(from_functions(model, scale))
  }

  stop("'model' must be a fitted gstat object, a named list of gstat variogram models, ",
       "a data frame of models (see availablemodels()) or a named list of functions", call. = F)
}

# FORMA DE LAS DISTANCIAS: Las funciones se evalúan sobre un vector y el
# resultado toma la forma de 'h' (el solver les pasa matrices de distancias)
keep_shape <- function(f){
  function(h){
    out <- if (length(h) == 0) numeric(0) else f(as.vector(h))
    dim(out) <- dim(h)
    out
  }
}

# MODELOS DE gstat: Se evalúan con gstat::variogramLine(), el mismo código con
# el que gstat los ajustó. La meseta es la covarianza en cero; gstat no la
# calcula para los modelos sin meseta (Lin con rango cero, Log y Pow).
from_gstat <- function(models){
  if (!requireNamespace('gstat', quietly = TRUE)) {
    stop("The package 'gstat' is needed to use gstat models", call. = F)
  }

  entries <- lapply(names(models), function(id) {
    m <- models[[id]]

    # La anisotropía cambia la distancia según la dirección, y el solver solo
    # trabaja con distancias
    if (any(m$anis1 != 1 | m$anis2 != 1)) {
      stop("Anisotropic variogram models are not supported yet by SpatFD (model of '", id,
           "'). Fit an isotropic model.", call. = F)
    }

    sill <- tryCatch(gstat::variogramLine(m, dist_vector = 0, covariance = TRUE)$gamma,
                     error = function(e) NA_real_)

    list(variogram  = keep_shape(function(h) gstat::variogramLine(m, dist_vector = h)$gamma),
         sill       = sill,
         structures = paste(as.character(m$model), collapse = ', '),
         source     = 'gstat',
         lmc        = data.frame(model = as.character(m$model), range = m$range,
                                 kappa = m$kappa, psill = m$psill))
  })

  stats::setNames(entries, names(models))
}

# MODELOS SIN MESETA DEL CATÁLOGO: Los mismos que en gstat. El lineal solo lo
# es con rango cero (o sin rango, que toma el valor por defecto, cero).
unbounded <- function(model, rango){
  model %in% c('logarithmic', 'power') | (model == 'linear' & (is.na(rango) | rango == 0))
}

# DATA FRAME DE MODELOS: Una fila por estructura. Cada estructura es el modelo
# del catálogo ('models') con sus parámetros fijos; los parámetros en NA quedan
# en el valor por defecto del modelo. El variograma de cada id suma sus
# estructuras.
from_table <- function(df){
  stopifnot(
    "The data frame of models needs the columns 'var', 'model' and 'sigma'" = all(c('var', 'model', 'sigma') %in% colnames(df)),
    "The data frame of models is empty" = nrow(df) > 0
  )

  df$model <- as.character(df$model)
  unknown  <- setdiff(df$model, names(models))
  if (length(unknown) > 0) {
    stop('Unknown variogram models: ', paste(unknown, collapse = ', '),
         '. See availablemodels(txt = TRUE)', call. = F)
  }

  params <- intersect(c('rango', 'sigma', 'alpha1'), colnames(df))
  ids    <- unique(as.character(df$var))

  entries <- lapply(ids, function(id) {
    rows <- df[df$var == id, , drop = F]

    structures <- lapply(seq_len(nrow(rows)), function(i) {
      model <- models[[rows$model[i]]]
      args  <- as.list(rows[i, params, drop = F])
      args  <- args[!vapply(args, is.na, logical(1))]

      component <- function(h) do.call(model, c(list(h), args))
      component(numeric(0))          # Parámetros inválidos fallan aquí y no al predecir
      component
    })

    rango   <- if ('rango' %in% colnames(rows)) rows$rango else rep(NA, nrow(rows))
    bounded <- !any(unbounded(rows$model, rango))

    list(variogram  = keep_shape(function(h) Reduce(`+`, lapply(structures, function(component) component(h)))),
         sill       = if (bounded) sum(rows$sigma) else NA_real_,
         structures = paste(rows$model, collapse = ', '),
         source     = 'data frame',
         lmc        = data.frame(model = rows$model, range = rango,
                                 kappa = if ('alpha1' %in% colnames(rows)) rows$alpha1 else NA,
                                 psill = rows$sigma))
  })

  stats::setNames(entries, ids)
}

# FUNCIONES DEL USUARIO: Cada una es un variograma de la distancia, su único
# argumento. La meseta se calcula una sola vez, al crear el modelo, ver
# function_sill().
from_functions <- function(functions, scale){
  entries <- lapply(names(functions), function(id) {
    f <- functions[[id]]

    if (length(formals(f)) != 1) {
      stop("The function of '", id, "' must have a single argument, the distance", call. = F)
    }

    test <- tryCatch(f(c(0, 1, 2)), error = function(e) NULL)
    if (!is.numeric(test) || length(test) != 3 || anyNA(test) || test[1] != 0) {
      stop("The function of '", id, "' must be a vectorised variogram: one number per ",
           'distance and 0 at distance 0', call. = F)
    }

    list(variogram  = keep_shape(f),
         sill       = function_sill(f, id, scale),
         structures = 'function',
         source     = 'function')
  })

  stats::setNames(entries, names(functions))
}

# MESETA DE UNA FUNCIÓN: Como gamma(h) = C(0) - C(h) y la covarianza se apaga
# (o, si oscila, promedia cero) a grandes distancias, la meseta C(0) es el
# valor de gamma lejos del origen. Se busca en cascada:
#   1. f(Inf)            -> finito: la meseta; infinito: modelo sin meseta
#   2. Límite numérico   -> si f se estabiliza a grandes distancias (wave, hole)
#   3. Promedio          -> si f oscila sin crecer (periódico): el promedio de f
#                           a grandes distancias. Si crece, no tiene meseta
# Las distancias de prueba van de 100 a 10^12 veces el tamaño del área de los
# sitios, en escala logarítmica: cada punto cae en una fase distinta de una
# oscilación, así que el promedio no depende del período.
function_sill <- function(f, id, scale){
  no_sill <- function(why) {
    message("The variogram of '", id, "' ", why, ': it is used as a model without a sill ',
            '(ordinary kriging only)')
    NA_real_
  }

  # 1. EN EL INFINITO
  value <- tryCatch(suppressWarnings(f(Inf)), error = function(e) NaN)
  if (is.numeric(value) && length(value) == 1 && !is.na(value)) {
    return(if (is.finite(value)) value else NA_real_)
  }

  x <- scale * 10^seq(2, 12, length.out = 1e5)
  y <- tryCatch(suppressWarnings(f(x)), error = function(e) NULL)
  if (!is.numeric(y) || length(y) != length(x) || any(!is.finite(y))) {
    return(no_sill('cannot be evaluated at large distances'))
  }

  # 2. LÍMITE NUMÉRICO: Los valores más lejanos ya no cambian
  far <- y[x >= scale * 1e10]
  if (diff(range(far)) <= 1e-6 * max(abs(far), .Machine$double.eps)) {
    message("The sill of '", id, "' is its limit at large distances")
    return(mean(far))
  }

  # 3. OSCILACIÓN: Acotada si la parte lejana no supera a la cercana
  near <- y[x < scale * 1e7]
  if (max(abs(y[x >= scale * 1e7])) > 1.5 * max(abs(near))) {
    return(no_sill('grows without a limit'))
  }

  message("The variogram of '", id, "' oscillates: its sill is its mean value at large distances")
  mean(y)
}

# ============================== IDS DEL MODELO ================================
# SCORES DE LOS OBJETOS: Una fila por score disponible
score_table <- function(objects){
  rows <- lapply(names(objects), function(name) {
    K <- ncol(get_pca(objects[[name]])$scores)
    data.frame(id = paste0(name, '.PC', seq_len(K)), variable = name, pc = seq_len(K))
  })
  do.call(rbind, rows)
}

# EMPAREJAMIENTO: Cada id del modelo debe ser un score (variograma directo) o
# una pareja de scores en cualquier orden (variograma cruzado). Se generan
# todos los ids posibles y se buscan coincidencias exactas: así funcionan
# nombres con puntos, como 'PM2.5'.
match_ids <- function(ids, objects){
  scores <- score_table(objects)

  pairs <- expand.grid(first = scores$id, second = scores$id, stringsAsFactors = F)
  pairs <- pairs[pairs$first != pairs$second, ]
  known <- rbind(data.frame(key = scores$id, first = scores$id, second = scores$id),
                 data.frame(key = paste(pairs$first, pairs$second, sep = '.'), pairs))

  # Nombres de variables muy particulares podrían dar dos lecturas del mismo id
  ambiguous <- intersect(ids, known$key[duplicated(known$key)])
  if (length(ambiguous) > 0) {
    stop('Ambiguous variogram ids: ', paste(ambiguous, collapse = ', '),
         '. Rename the variables with the argument name of spatfd()', call. = F)
  }

  unknown <- setdiff(ids, known$key)
  if (length(unknown) > 0) {
    stop('The variogram ids do not match the scores of the objects: ', paste(unknown, collapse = ', '),
         '\nThe scores of the objects are: ', paste(scores$id, collapse = ', '), call. = F)
  }

  matched <- known[match(ids, known$key), ]
  data.frame(id     = ids,
             type   = ifelse(matched$first == matched$second, 'direct', 'cross'),
             first  = matched$first,
             second = matched$second,
             row.names = NULL)
}

# NÚMERO DE PC DE CADA VARIABLE: Los variogramas directos de cada variable
# deben ser los de sus primeras K componentes, sin saltos. Los cruzados solo
# pueden usar esas componentes.
deduce_k <- function(ids, objects){
  scores <- score_table(objects)
  direct <- scores[scores$id %in% ids$id[ids$type == 'direct'], ]

  K <- vapply(names(objects), function(name) {
    pcs <- direct$pc[direct$variable == name]

    if (length(pcs) == 0) {
      stop("The model has no variogram for the variable '", name, "'. Every variable given needs ",
           "the variograms of its first components (", name, ".PC1, ...)", call. = F)
    }
    if (!setequal(pcs, seq_len(max(pcs)))) {
      stop("The variograms of '", name, "' must be for its first components, without gaps. Found: ",
           paste0('PC', sort(pcs), collapse = ', '), call. = F)
    }
    max(pcs)
  }, integer(1))

  used  <- c(ids$first, ids$second)[ids$type == 'cross']
  extra <- setdiff(used, direct$id)
  if (length(extra) > 0) {
    stop('Cross-variograms with scores that have no direct variogram: ', paste(unique(extra), collapse = ', '),
         call. = F)
  }

  K
}

# ============================ COPIAS TRUNCADAS ================================
# COPIA DE UNA VARIABLE CON SUS PRIMERAS K COMPONENTES: Es la que se guarda en
# el modelo para predecir y validar. Las curvas, los datos y los tiempos se
# conservan; de la PCA solo quedan las K componentes del modelo.
truncate_spatfd <- function(sfd, K){
  name <- get_name(sfd)
  pca  <- get_pca(sfd)
  keep <- seq_len(K)

  pca$harmonics <- pca$harmonics[keep]
  pca$scores    <- pca$scores[, keep, drop = F]
  pca$varprop   <- pca$varprop[keep]

  internal     <- attr(sfd, 'internal')
  internal$pca <- pca

  coords   <- get_coords(sfd)
  scores   <- sfd$scores[, c(colnames(coords), paste0(name, '.PC', keep)), drop = F]
  variance <- sfd$variance
  variance[['no. PC']]  <- K
  variance[['variance']] <- sum(pca$varprop)

  new_spatfd(scores, variance, sfd$missing, internal, name)
}

# ===================== VALIDEZ DE UN MODELO DE CORREGIONALIZACIÓN =============
# Cuando los variogramas de un grupo de scores conectados por cruzados usan las
# mismas estructuras (mismo modelo y rango), el modelo es un modelo lineal de
# corregionalización: es válido si, en cada estructura, la matriz de mesetas
# parciales de los scores es semidefinida positiva. Es el mismo criterio de
# gstat y vale para cualquier conjunto de sitios. Los modelos que no tienen esa
# forma (por ejemplo, funciones) se revisan solo en los sitios observados, al
# armar el sistema, ver solver_UTILS.R.
check_coregionalization <- function(compiled, ids, data){
  fields <- unlist(model_fields(data), use.names = F)
  matrix <- score_variograms(fields, ids, lapply(compiled, `[[`, 'variogram'))

  for (group in solver_groups(matrix)) {
    if (length(group) == 1) next

    members <- ids$id[ids$first %in% group & ids$second %in% group]
    tables  <- lapply(compiled[members], `[[`, 'lmc')
    if (any(vapply(tables, is.null, logical(1)))) next

    # Las estructuras de un variograma, como texto: modelo, rango y forma
    signature <- function(table) sort(paste(table$model, table$range, table$kappa))
    if (length(unique(lapply(tables, signature))) > 1) next

    for (structure in signature(tables[[1]])) {
      B <- matrix(0, length(group), length(group), dimnames = list(group, group))
      for (id in members) {
        row   <- ids[ids$id == id, ]
        table <- tables[[id]]
        B[row$first, row$second] <- B[row$second, row$first] <-
          table$psill[paste(table$model, table$range, table$kappa) == structure]
      }

      values <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
      if (min(values) < -1e-8 * max(abs(values))) stop('The model is not valid', call. = F)
    }
  }

  invisible(TRUE)
}

# ========================= SISTEMAS DE KRIGING ================================
# Al crear el modelo se arman, validan y factorizan los sistemas de los dos
# métodos de predicción, ver solver_UTILS.R. Así predict() solo resuelve.

# SCORES DE LAS COPIAS: Ids de los scores de cada variable, en orden
model_fields <- function(data){
  lapply(data, function(sfd) paste0(get_name(sfd), '.PC', seq_len(ncol(get_pca(sfd)$scores))))
}

# MATRIZ DE VARIOGRAMAS DE LOS SCORES: La que necesita el solver, con NULL en
# las parejas de scores sin variograma
score_variograms <- function(fields, ids, variograms){
  q   <- length(fields)
  out <- matrix(vector('list', q * q), q, q, dimnames = list(fields, fields))
  for (r in seq_len(nrow(ids))) {
    out[[ids$first[r], ids$second[r]]] <- variograms[[ids$id[r]]]
    out[[ids$second[r], ids$first[r]]] <- variograms[[ids$id[r]]]
  }
  out
}

# MÉTODO SCORES: Cada score es un campo, con los sitios de su variable
scores_systems <- function(data, ids, variograms){
  fields <- model_fields(data)
  sites  <- unlist(lapply(names(fields), function(v) {
    stats::setNames(rep(list(get_coords(data[[v]])), length(fields[[v]])), fields[[v]])
  }), recursive = F)

  matrix <- score_variograms(unlist(fields, use.names = F), ids, variograms)
  lapply(solver_groups(matrix), solver_system, sites = sites, variograms = matrix)
}

# MÉTODO LAMBDA: Cada variable es un campo, su curva. El variograma entre las
# curvas de las variables p y q se arma con los de sus scores,
#
#     gamma_pq(h) = suma_k suma_l  G[k, l] * gamma_{pk, ql}(h)
#
# donde G[k, l] es el producto interno entre las funciones propias k de p y l
# de q sobre el intervalo de tiempo común. Dentro de una variable observada en
# todo ese intervalo G es la identidad, porque sus funciones propias son
# ortonormales. Si las variables no comparten tiempos no hay sistema: el
# modelo sirve para el método 'scores', pero no para 'lambda'.
lambda_systems <- function(data, ids, variograms){
  fields <- model_fields(data)
  scores <- score_variograms(unlist(fields, use.names = F), ids, variograms)

  ranges <- vapply(data, function(sfd) get_fd(sfd)$basis$rangeval, numeric(2))
  domain <- c(max(ranges[1, ]), min(ranges[2, ]))
  if (domain[1] >= domain[2]) return(list(domain = NULL, systems = NULL))
  own <- ranges[1, ] == domain[1] & ranges[2, ] == domain[2]

  curve_variogram <- function(p, q){
    G <- if (p == q && own[[p]]) diag(length(fields[[p]])) else
      inner_product(get_pca(data[[p]])$harmonics, get_pca(data[[q]])$harmonics, domain)

    terms <- list()
    for (k in seq_along(fields[[p]])) {
      for (l in seq_along(fields[[q]])) {
        gamma <- scores[[fields[[p]][k], fields[[q]][l]]]
        if (!is.null(gamma) && G[k, l] != 0) terms <- c(terms, list(list(G = G[k, l], gamma = gamma)))
      }
    }
    if (length(terms) == 0) return(NULL)                 # Variables sin variograma cruzado
    function(h) Reduce(`+`, lapply(terms, function(term) term$G * term$gamma(h)))
  }

  variables <- names(data)
  matrix    <- matrix(vector('list', length(variables)^2), length(variables), length(variables),
                      dimnames = list(variables, variables))
  for (p in variables) for (q in variables) matrix[p, q] <- list(curve_variogram(p, q))

  sites <- lapply(data, get_coords)
  list(domain  = domain,
       systems = lapply(solver_groups(matrix), solver_system, sites = sites, variograms = matrix))
}

# ============================ ACCESORES DEL MODELO ============================
get_model_data <- function(model) attr(model, 'data')
get_compiled   <- function(model) attr(model, 'internal')

# VARIOGRAMA DE UNA PAREJA DE SCORES: En el orden que sea; NULL si el modelo no
# lo tiene (los scores se toman como independientes)
pair_variogram <- function(model, a, b){
  ids <- get_compiled(model)$ids
  row <- which((ids$first == a & ids$second == b) | (ids$first == b & ids$second == a))
  if (length(row) == 0) return(NULL)
  get_compiled(model)$variograms[[ids$id[row]]]
}

# ========================== GRÁFICO DE MODELOS ================================
# VARIABLES DE VARIOS MODELOS: Para el gráfico, cada variable que aparece en
# alguno de los modelos, con la copia que tiene más componentes. Si una variable
# con el mismo nombre tiene scores distintos en dos modelos, se avisa.
models_data <- function(models){
  data <- list()

  for (model in models) {
    for (copy in get_model_data(model)) {
      name <- get_name(copy)
      if (is.null(data[[name]])) {
        data[[name]] <- copy
        next
      }

      old <- get_pca(data[[name]])$scores
      new <- get_pca(copy)$scores
      k   <- seq_len(min(ncol(old), ncol(new)))
      if (!isTRUE(all.equal(old[, k, drop = F], new[, k, drop = F], check.attributes = F))) {
        warning("The variable '", name, "' has different scores in the models. The empirical ",
                'variograms use the first model where it appears', call. = F)
      } else if (ncol(new) > ncol(old)) {
        data[[name]] <- copy
      }
    }
  }

  data
}

# VARIOGRAMAS EMPÍRICOS: Se calculan con gstat::variogram() sobre los scores
# guardados en los modelos, para que sean los mismos puntos que el usuario vio
# al ajustar. Se devuelven con las parejas de scores de cada fila.
empirical_variograms <- function(data, ...){
  K <- vapply(data, function(sfd) ncol(get_pca(sfd)$scores), integer(1))
  g <- do.call(to_gstat, c(unname(data), list(pcs = K)))
  v <- as.data.frame(gstat::variogram(g, ...))

  # gstat nombra los cruzados 'a.b' con a antes que b en el objeto gstat
  fields <- names(g$data)
  pairs  <- expand.grid(first = fields, second = fields, stringsAsFactors = F)
  key    <- ifelse(pairs$first == pairs$second, pairs$first, paste(pairs$first, pairs$second, sep = '.'))
  found  <- match(as.character(v$id), key)

  data.frame(dist = v$dist, gamma = v$gamma, first = pairs$first[found], second = pairs$second[found])
}
