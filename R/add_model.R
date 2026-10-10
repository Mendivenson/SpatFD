#' Add a variogram model to a SpatFD object
#'
#' Registers a set of fitted variogram models inside a `SpatFD` object under a
#' name, ready to be used by [predict()][predict.SpatFD]. An object can hold
#' several models, so that they can be compared later; each one is compiled
#' only once, when it is added.
#'
#' @param sfd An object of class `SpatFD`, as returned by [spatfd()].
#' @param model The models to add. Either fitted models that [toSpatFD()] can
#'   translate (a `gstat` object or a named list of `variogramModel` objects),
#'   or a data frame of models as returned by [toSpatFD()].
#'
#'   The data frame must have the columns `var` (id of the variogram), `model`
#'   (a name in [models]) and `sigma`; every column other than `var` and
#'   `model` is taken as a model parameter (`rango`, `sigma`, `alpha1`).
#'   Parameters set to `NA` are left at the default of the model, and
#'   parameters outside the domain of their model raise an error (see
#'   [models]). The ids in `var` must be scores of `sfd` (`<variable>.PC1`,
#'   ...) for direct variograms, or two of them joined by a dot for
#'   cross-variograms.
#' @param name Name given to the model. By default the models are named
#'   sequentially: `"model 1"`, `"model 2"`, ... If the name is already in
#'   use, that model is replaced.
#'
#' @details
#' The model is added to `sfd` **in place**: there is no need to assign the
#' result. The compiled models live in the attribute `models`, an environment
#' like the other attributes of the object (see [spatfd()]), so copies of the
#' object made with `<-` share the same models.
#'
#' `sfd$models` shows a summary table of the registered models, with one row
#' per model and the columns
#' \describe{
#'   \item{`variables`}{Functional variables the model covers.}
#'   \item{`structures`}{Basic variogram models it is made of.}
#'   \item{`no. scores`}{Number of scores with a direct variogram.}
#'   \item{`cross`}{Whether it includes cross-variograms.}
#'   \item{`sill`}{Sum of the sills of the direct variograms.}
#' }
#'
#' Each entry of `attr(sfd, "models")` is a list with the table of the model
#' (`table`), the compiled variograms (`variograms`, a named list of functions
#' of the distance `h`) and their sills (`sills`).
#'
#' @return `sfd`, invisibly.
#'
#' @section Theory:
#' A variogram is usually modelled as a nested sum of basic structures, each
#' with its own partial sill \eqn{\sigma_j} and range \eqn{a_j}, often
#' including a nugget effect:
#' \deqn{\gamma(h) = \sum_{j=1}^{J} \gamma_j(h; \sigma_j, a_j)}
#' Each row of the table of models is one structure \eqn{\gamma_j}, and the compiled function
#' is their sum \eqn{\gamma(h)} (Chilès and Delfiner, 2012).
#'
#' @references
#' Chilès, J.-P. and Delfiner, P. (2012). *Geostatistics: Modeling Spatial
#' Uncertainty* (2nd ed.). John Wiley & Sons. \doi{10.1002/9781118136188}
#'
#' @seealso [toSpatFD()], [models], [availablemodels()], [predict.SpatFD()]
#'
#' @examples
#' set.seed(1)
#' t <- 1:50
#' dat <- sapply(1:10, function(i) rnorm(1, 5) * sin(2 * pi * t / 50) + rnorm(50, sd = 0.1))
#' colnames(dat) <- paste0("s", 1:10)
#' coords <- cbind(x = runif(10), y = runif(10))
#' rownames(coords) <- colnames(dat)
#' sfd <- spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "sim")
#'
#' df <- data.frame(
#'   var    = c("sim.PC1", "sim.PC1", "sim.PC2"),
#'   model  = c("nugget", "exponential", "spherical"),
#'   rango  = c(NA, 0.5, 0.8),
#'   sigma  = c(2, 30, 0.01),
#'   alpha1 = c(NA, NA, NA)
#' )
#'
#' add_model(sfd, df)
#' add_model(sfd, df[df$model != "nugget", ], name = "no nugget")
#'
#' # Fitted models of gstat are translated on the way in
#' fitted <- list(
#'   "sim.PC1" = gstat::vgm(psill = 30, model = "Sph", range = 0.5),
#'   "sim.PC2" = gstat::vgm(psill = 0.01, model = "Sph", range = 0.5)
#' )
#' add_model(sfd, fitted, name = "spherical")
#'
#' sfd$models
#'
#' vg <- attr(sfd, "models")[["model 1"]]$variograms
#' vg$sim.PC1(c(0, 0.1, 0.5, 2))
#'
#' @export
add_model <- function(sfd, model, name = NULL){

  stopifnot(
    'function only available for SpatFD objects' = inherits(sfd, 'SpatFD'),
    "'name' must be a single string or null" =
      is.null(name) || (is.character(name) && length(name) == 1 && !is.na(name) && nzchar(name))
  )

  # Lo que no venga ya en el esquema tabular de SpatFD se traduce primero
  is_table <- is.data.frame(model) && all(c('var', 'model', 'sigma') %in% colnames(model))
  df <- if (is_table) model else toSpatFD(model)

  stopifnot('There are no models to add' = nrow(df) > 0)

  # Los ids deben ser scores del objeto o parejas de scores (variogramas cruzados)
  ids     <- names(score_ids(sfd))
  unknown <- setdiff(df$var, c(ids, outer(ids, ids, paste, sep = '.')))
  if (length(unknown) > 0) {
    stop('The variogram ids do not match the scores of the object: ',
         paste(unknown, collapse = ', '),
         '\nThe scores of the object are: ', paste(ids, collapse = ', '), call. = FALSE)
  }

  unknown <- setdiff(df$model, names(models))
  if (length(unknown) > 0) {
    stop('Unknown variogram models: ', paste(unknown, collapse = ', '),
         '. See availablemodels(txt = TRUE)', call. = FALSE)
  }

  # REGISTRO: El modelo se guarda en el entorno del atributo 'models', por lo
  # que el objeto se modifica en el sitio
  registry <- attr(sfd, 'models')
  index    <- vapply(as.list(registry), function(m) m$index, numeric(1))

  if (is.null(name)) {
    name <- next_model_name(registry)
  } else if (name %in% names(index)) {
    message("Replacing the model '", name, "'")
  }

  registry[[name]] <- list(
    index      = if (name %in% names(index)) index[[name]] else max(0, index) + 1,
    table      = df,
    variograms = compile_model(df),
    sills      = vapply(split(df$sigma, df$var), sum, numeric(1))
  )

  invisible(sfd)
}

# .............................. FUNCIONES INTERNAS ................................................

# NOMBRES SECUENCIALES: model 1, model 2, ... saltando los nombres ya usados
next_model_name <- function(registry){
  k <- length(ls(registry)) + 1
  while (exists(paste('model', k), envir = registry, inherits = FALSE)) k <- k + 1
  paste('model', k)
}

# SCORES DEL OBJETO: Vector con la variable a la que pertenece cada score,
# nombrado con el id del score (<variable>.PC1, <variable>.PC2, ...)
score_ids <- function(sfd){
  variables <- var_names(sfd)
  ids <- lapply(variables, function(v) {
    paste0(v, '.PC', seq_len(ncol(attr(sfd, 'pca')[[v]]$scores)))
  })
  stats::setNames(rep(variables, lengths(ids)), unlist(ids))
}

# TABLA RESUMEN DE LOS MODELOS REGISTRADOS: Es lo que el usuario ve en sfd$models
model_table <- function(sfd){
  ids     <- score_ids(sfd)
  entries <- as.list(attr(sfd, 'models'))
  entries <- entries[order(vapply(entries, function(m) m$index, numeric(1)))]

  rows <- lapply(entries, function(m) {
    direct <- intersect(names(ids), m$table$var)
    data.frame('variables'  = paste(unique(ids[direct]), collapse = ', '),
               'structures' = paste(unique(m$table$model), collapse = ', '),
               'no. scores' = length(direct),
               'cross'      = !all(m$table$var %in% names(ids)),
               'sill'       = sum(m$sills[direct]),
               check.names = F)
  })

  if (length(rows) == 0) {
    return(data.frame('variables' = character(0), 'structures' = character(0),
                      'no. scores' = integer(0), 'cross' = logical(0),
                      'sill' = numeric(0), check.names = F))
  }

  rows <- do.call(rbind, rows)
  rownames(rows) <- names(entries)
  rows
}

# COMPILACIÓN: Una función de la distancia por cada id de variograma en df$var
compile_model <- function(df){

  # Función para fijar los argumentos de una función dada
  fix_args <- function(f, args_to_fix) {
    new_f <- f
    body(new_f) <- do.call(substitute, list(body(f), args_to_fix))
    formals(new_f) <- formals(f)[!names(formals(f)) %in% names(args_to_fix)]
    environment(new_f) <- globalenv()
    new_f
  }

  # Función para forzar la suma de funciones como una nueva función
  combine_models <- function(componentes) {
    llamadas <- lapply(componentes, function(fn) as.call(list(fn, quote(h))))
    expr_suma <- Reduce(function(a, b) bquote(.(a) + .(b)), llamadas)

    f <- function(h) NULL
    body(f) <- expr_suma
    environment(f) <- globalenv()
    f
  }

  all_models = list()
  for (v in unique(df$var)){
    model <- df[df$var == v,]

    componentes <- lapply(seq_len(nrow(model)), function(i) {
      x <- as.list(model[i, ])
      m <- x$model
      params <- x[setdiff(names(x), c('model', 'var'))]
      params <- params[!vapply(params, is.na, logical(1))]  # quitar NA (params que ese modelo no usa)
      componente <- fix_args(models[[m]], params)
      componente(numeric(0))  # Parámetros inválidos fallan aquí y no al evaluar
      componente
    })

    all_models[[v]] = combine_models(componentes)
  }
  all_models
}
