# ====================== FUNCIONES INTERNAS DE predict() =======================
# Solo las usa predict() (y sus métodos). El sistema de kriging está en
# solver_UTILS.R; aquí se le dan los valores observados y se pasa de lo que
# predice el solver a curvas.
#
# Los dos métodos devuelven, para cada variable pedida:
#   fd   -> La curva predicha en cada sitio nuevo, como fd (para los gráficos)
#   pred -> Las curvas evaluadas en los tiempos observados (tiempos x sitios)
#   var  -> La varianza de predicción en cada sitio nuevo

# ============================ SITIOS NUEVOS ===================================
# Matriz con las coordenadas de los sitios nuevos: columnas por nombre si los
# tiene, por posición si no. Los sitios sin nombre se llaman new1, new2, ...
check_newcoords <- function(newcoords, coord_names){
  stopifnot(
    "'newcoords' must be a numeric matrix or data frame without missing values" =
      (is.matrix(newcoords) || is.data.frame(newcoords)) && nrow(newcoords) > 0 &&
      is.numeric(as.matrix(newcoords)) && !anyNA(newcoords)
  )

  newcoords <- as.matrix(newcoords)
  if (all(coord_names %in% colnames(newcoords))) {
    newcoords <- newcoords[, coord_names, drop = F]
  } else if (is.null(colnames(newcoords)) && ncol(newcoords) == length(coord_names)) {
    colnames(newcoords) <- coord_names
  } else {
    stop("'newcoords' must have the coordinates of the model: ", paste(coord_names, collapse = ', '), call. = F)
  }

  if (is.null(rownames(newcoords))) rownames(newcoords) <- paste0('new', seq_len(nrow(newcoords)))
  newcoords
}

# GRILLA SOBRE UN MAPA: Cerca de grid_size puntos en una grilla regular dentro del mapa
map_grid <- function(map, grid_size){
  if (!requireNamespace('sf', quietly = TRUE)) {
    stop("The package 'sf' is needed to predict over a map", call. = F)
  }
  stopifnot("'grid_size' must be a single positive number" = is.numeric(grid_size) && length(grid_size) == 1 && grid_size > 0)

  area <- sf::st_union(sf::st_geometry(map))
  grid <- sf::st_coordinates(sf::st_sample(area, size = grid_size, type = 'regular'))[, 1:2, drop = F]
  dimnames(grid) <- NULL
  grid
}

# TIEMPOS COMO NOMBRES DE FILA: Los del usuario, como texto
time_labels <- function(time){
  if (is.numeric(time)) as.character(time) else format(time)
}

# ============================= MÉTODO SCORES ==================================
# KRIGING DE LOS SCORES: Predicción y varianza de cada score en cada sitio
# nuevo (sitios nuevos x scores). Cada grupo de scores conectados se resuelve
# con su sistema, ya factorizado al crear el modelo.
krige_scores <- function(object, newcoords){
  data   <- get_model_data(object)
  fields <- unlist(model_fields(data), use.names = F)

  # Valores observados de cada score: una fila, una columna por sitio
  values <- unlist(lapply(data, function(sfd) {
    scores <- get_pca(sfd)$scores
    lapply(seq_len(ncol(scores)), function(k) matrix(scores[, k], nrow = 1))
  }), recursive = F)
  names(values) <- fields

  pred <- variance <- matrix(NA_real_, nrow(newcoords), length(fields),
                             dimnames = list(rownames(newcoords), fields))

  for (system in get_compiled(object)$scores) {
    solution <- solver_solve(system, newcoords)
    for (field in system$fields) {
      pred[, field]     <- solver_combine(system, solution[[field]], values)
      variance[, field] <- solution[[field]]$variance
    }
  }

  list(pred = pred, variance = variance)
}

# CURVAS: Media más los scores predichos por las funciones propias. Como las
# funciones propias son ortonormales, la varianza de la curva es la suma de las
# varianzas de sus scores.
predict_scores <- function(object, newcoords, variable){
  kriging <- krige_scores(object, newcoords)
  data    <- get_model_data(object)
  fields  <- model_fields(data)

  curves <- lapply(variable, function(v) {
    pca <- get_pca(data[[v]])

    coefs <- pca$harmonics$coefs %*% t(kriging$pred[, fields[[v]], drop = F]) + as.vector(pca$meanfd$coefs)
    colnames(coefs) <- rownames(newcoords)
    fd <- fda::fd(coefs, pca$harmonics$basis)

    pred <- fda::eval.fd(get_internal(data[[v]], 't'), fd)
    dimnames(pred) <- list(time_labels(get_internal(data[[v]], 'time')), rownames(newcoords))

    list(fd = fd, pred = pred, var = rowSums(kriging$variance[, fields[[v]], drop = F]))
  })

  stats::setNames(curves, variable)
}

# ============================= MÉTODO LAMBDA ==================================
# La curva de cada variable es la suma de las curvas observadas por sus pesos.
# Con varias variables, en los tiempos que comparten (el intervalo común).
predict_lambda <- function(object, newcoords, variable){
  lambda <- get_compiled(object)$lambda
  if (is.null(lambda$domain)) {
    stop("The variables of the model have no time in common, so the method 'lambda' cannot ",
         'combine them', call. = F)
  }

  data   <- get_model_data(object)
  domain <- lambda$domain

  curves <- lapply(variable, function(r) {
    system   <- Find(function(s) r %in% s$fields, lambda$systems)
    solution <- solver_solve(system, newcoords)[[r]]

    # Tiempos de la variable dentro del intervalo común
    t    <- get_internal(data[[r]], 't')
    keep <- t >= domain[1] & t <= domain[2]
    t    <- t[keep]

    # Curvas observadas de cada variable del grupo, evaluadas en esos tiempos
    values <- lapply(stats::setNames(nm = system$fields), function(p) fda::eval.fd(t, get_fd(data[[p]])))

    pred <- solver_combine(system, solution, values)
    dimnames(pred) <- list(time_labels(get_internal(data[[r]], 'time')[keep]), rownames(newcoords))

    # Versión funcional, solo para graficar: en la base de la variable si se
    # predijo en todo su intervalo, en B-splines sobre el intervalo común si no
    own   <- all(get_fd(data[[r]])$basis$rangeval == domain)
    basis <- if (own) get_fd(data[[r]])$basis else
      fda::create.bspline.basis(domain, get_internal(data[[r]], 'args')$nbasis)
    fd <- fda::smooth.basis(t, pred, basis)$fd
    colnames(fd$coefs) <- rownames(newcoords)

    list(fd = fd, pred = pred, var = solution$variance)
  })

  if (!all(vapply(data, function(sfd) all(get_fd(sfd)$basis$rangeval == domain), logical(1)))) {
    message("Method 'lambda': the curves are predicted over the time interval common to the variables")
  }

  stats::setNames(curves, variable)
}
