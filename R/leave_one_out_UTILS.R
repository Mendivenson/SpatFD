# ===================== FUNCIONES INTERNAS DE leave_one_out() ==================
# Solo las usa leave_one_out() (y su método print).

# ETIQUETAS DE LOS MODELOS: El nombre del argumento si lo tiene, o la expresión
# con la que se dio el modelo
model_labels <- function(models, exprs){
  given  <- if (is.null(names(models))) rep('', length(models)) else names(models)
  labels <- ifelse(given == '', vapply(exprs, function(e) paste(deparse(e), collapse = ''), character(1)), given)
  make.unique(labels)
}

# VARIABLE SIN UN SITIO: Repite spatfd() sin el sitio, con los mismos datos,
# tiempos y argumentos del suavizamiento, y con las mismas K componentes del
# modelo. Así la curva que se deja fuera no participa ni en la media ni en las
# componentes con que se predice.
refit_without <- function(copy, site){
  args  <- get_internal(copy, 'args')
  data  <- get_internal(copy, 'data')
  keep  <- colnames(data) != site

  dat <- data.frame(.time = get_internal(copy, 'time'), data[, keep, drop = F], check.names = F)
  sfd <- suppressWarnings(suppressMessages(
    spatfd(dat, get_coords(copy)[keep, , drop = F], time = '.time',
           basis = args$basis, nbasis = args$nbasis, lambda = args$lambda,
           nharm = ncol(get_pca(copy)$scores), name = get_name(copy))
  ))

  align_signs(sfd, copy)
}

# SIGNOS DE LAS COMPONENTES: El signo de una componente principal es
# arbitrario: al repetir el análisis sin un sitio, una componente puede salir
# con el signo contrario, y entonces sus variogramas cruzados cambiarían de
# sentido. Cada componente toma el signo de la del análisis completo: si su
# producto interno con ella es negativo, se cambia el signo de la componente y
# de sus scores.
align_signs <- function(sfd, reference){
  new <- get_pca(sfd)
  old <- get_pca(reference)

  range <- c(max(new$harmonics$basis$rangeval[1], old$harmonics$basis$rangeval[1]),
             min(new$harmonics$basis$rangeval[2], old$harmonics$basis$rangeval[2]))
  flip  <- diag(inner_product(new$harmonics, old$harmonics, range)) < 0

  if (any(flip)) {
    name <- get_name(sfd)
    new$harmonics$coefs[, flip] <- -new$harmonics$coefs[, flip]
    new$scores[, flip]          <- -new$scores[, flip]

    internal     <- attr(sfd, 'internal')
    internal$pca <- new
    attr(sfd, 'internal') <- internal

    columns <- paste0(name, '.PC', which(flip))
    sfd$scores[, columns] <- -sfd$scores[, columns]
  }

  sfd
}

# CURVA DE ERROR DE UNA PREDICCIÓN: e(t) = chi(t) - pred(t), la diferencia
# entre la curva observada (suavizada) y la predicha, en los tiempos de la
# predicción (los de los datos). Se resume con la raíz del error cuadrático
# medio en el tiempo,
#
#     RMSE = sqrt( (1 / |T|) * integral_T e(t)^2 dt )
#
# que se calcula como sqrt(mean(e^2)): el promedio sobre los tiempos de los
# datos (con tiempos equiespaciados es exactamente esa integral). Queda en las
# unidades de los datos.
error_curve <- function(observed, predicted, t){
  as.vector(fda::eval.fd(t, observed)) - as.vector(predicted)
}

# TABLA RESUMEN DE UNA VARIABLE: Media, desviación estándar y máximo de los
# errores de los sitios, por modelo y método, de menor a mayor error medio
summary_table <- function(errors){
  keys <- unique(errors[, c('model', 'method')])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    e <- errors$error[errors$model == keys$model[i] & errors$method == keys$method[i]]
    data.frame(model = keys$model[i], method = keys$method[i], mean = mean(e), sd = stats::sd(e), max = max(e))
  })
  table <- do.call(rbind, rows)
  table <- table[order(table$mean), ]
  rownames(table) <- NULL
  table
}
