# ============================= FUNCIONES INTERNAS =============================
# Todo lo funcional del objeto (curvas suavizadas, componentes principales,
# datos originales y modelos compilados) vive en un único atributo oculto,
# 'internal', que es una lista normal:
#
#   coords    -> Coordenadas de todos los sitios del objeto
#   variables -> Una entrada por variable con: data, fd, pca y args
#   models    -> Una entrada por modelo con: table, variograms y sills
#
# El resto del paquete accede a él solo a través de estas funciones, así que
# el usuario nunca tiene que manipular objetos del paquete fda.

# CONSTRUCTOR: Arma el objeto a partir de sus partes. La tabla de modelos que
# ve el usuario se calcula aquí, así que nunca queda desactualizada.
new_spatfd <- function(scores, variance, coords, variables, models = list()){
  sfd <- structure(list(scores = scores, variance = variance, models = NULL),
                   class = 'SpatFD',
                   internal = list(coords = coords, variables = variables, models = models))
  sfd$models <- model_table(sfd)
  sfd
}

# ACCESORES
get_coords   <- function(sfd) attr(sfd, 'internal')$coords
get_variable <- function(sfd, variable) attr(sfd, 'internal')$variables[[variable]]
get_fd       <- function(sfd, variable) get_variable(sfd, variable)$fd
get_pca      <- function(sfd, variable) get_variable(sfd, variable)$pca
get_models   <- function(sfd) attr(sfd, 'internal')$models

# MODIFICADORES: Devuelven el objeto modificado, nunca lo cambian en el sitio
set_pca <- function(sfd, variable, pca){
  internal <- attr(sfd, 'internal')
  internal$variables[[variable]]$pca <- pca
  attr(sfd, 'internal') <- internal
  sfd
}

set_models <- function(sfd, models){
  internal <- attr(sfd, 'internal')
  internal$models <- models
  attr(sfd, 'internal') <- internal
  sfd$models <- model_table(sfd)
  sfd
}

# SCORES DEL OBJETO: Vector con la variable a la que pertenece cada score,
# nombrado con el id del score (<variable>.PC1, <variable>.PC2, ...)
score_ids <- function(sfd){
  variables <- var_names(sfd)
  ids <- lapply(variables, function(v) {
    paste0(v, '.PC', seq_len(ncol(get_pca(sfd, v)$scores)))
  })
  stats::setNames(rep(variables, lengths(ids)), unlist(ids))
}

# NOMBRES SECUENCIALES: model 1, model 2, ... saltando los nombres ya usados
next_model_name <- function(used){
  k <- length(used) + 1
  while (paste('model', k) %in% used) k <- k + 1
  paste('model', k)
}

# TABLA RESUMEN DE LOS MODELOS REGISTRADOS: Es lo que el usuario ve en sfd$models
model_table <- function(sfd){
  ids     <- score_ids(sfd)
  entries <- get_models(sfd)

  if (length(entries) == 0) {
    return(data.frame('variables' = character(0), 'structures' = character(0),
                      'no. scores' = integer(0), 'cross' = logical(0),
                      'sill' = numeric(0), check.names = F))
  }

  rows <- lapply(entries, function(m) {
    direct <- intersect(names(ids), m$table$var)
    data.frame('variables'  = paste(unique(ids[direct]), collapse = ', '),
               'structures' = paste(unique(m$table$model), collapse = ', '),
               'no. scores' = length(direct),
               'cross'      = !all(m$table$var %in% names(ids)),
               'sill'       = sum(m$sills[direct]),
               check.names = F)
  })

  rows <- do.call(rbind, rows)
  rownames(rows) <- names(entries)
  rows
}

# TABLA DE MODELOS PARA LOS MENSAJES DE ERROR
models_message <- function(sfd){
  paste(utils::capture.output(print(sfd$models)), collapse = '\n')
}

# CURVAS EVALUADAS: Matriz con una fila por tiempo y una columna por curva. Es
# la forma en que el usuario recibe las curvas, igual a la de los datos que
# entregó a spatfd().
evaluate_curves <- function(fd, t){
  values <- fda::eval.fd(t, fd)
  dimnames(values) <- list(t, colnames(fd$coefs))
  values
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
