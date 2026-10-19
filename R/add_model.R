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
#' Like any other modelling function in R, `add_model()` returns the updated
#' object, so the result must be assigned: `sfd <- add_model(sfd, model)`.
#'
#' `sfd$models` is the summary table of the registered models, with one row
#' per model and the columns
#' \describe{
#'   \item{`variables`}{Functional variables the model covers.}
#'   \item{`structures`}{Basic variogram models it is made of.}
#'   \item{`no. scores`}{Number of scores with a direct variogram.}
#'   \item{`cross`}{Whether it includes cross-variograms.}
#'   \item{`sill`}{Sum of the sills of the direct variograms.}
#' }
#'
#' @return `sfd` with the model added.
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
#' sfd <- add_model(sfd, df)
#' sfd <- add_model(sfd, df[df$model != "nugget", ], name = "no nugget")
#'
#' # Fitted models of gstat are translated on the way in
#' fitted <- list(
#'   "sim.PC1" = gstat::vgm(psill = 30, model = "Sph", range = 0.5),
#'   "sim.PC2" = gstat::vgm(psill = 0.01, model = "Sph", range = 0.5)
#' )
#' sfd <- add_model(sfd, fitted, name = "spherical")
#'
#' sfd$models
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

  # REGISTRO: Un nombre repetido reemplaza el modelo sin cambiar su posición
  registry <- get_models(sfd)
  if (is.null(name)) {
    name <- next_model_name(names(registry))
  } else if (name %in% names(registry)) {
    message("Replacing the model '", name, "'")
  }

  registry[[name]] <- list(
    table      = df,
    variograms = compile_model(df),
    sills      = vapply(split(df$sigma, df$var), sum, numeric(1))
  )

  set_models(sfd, registry)
}

# .............................. FUNCIONES INTERNAS ................................................

# COMPILACIÓN: Una función de la distancia por cada id de variograma en
# df$var, que suma las estructuras de ese variograma. Cada estructura es el
# modelo de 'models' con sus parámetros fijos; los parámetros en NA quedan en
# el valor por defecto del modelo.
compile_model <- function(df){
  params <- setdiff(colnames(df), c('var', 'model'))
  ids    <- unique(df$var)

  variograms <- lapply(ids, function(v) {
    rows <- df[df$var == v, , drop = FALSE]

    structures <- lapply(seq_len(nrow(rows)), function(i) {
      model <- models[[as.character(rows$model[i])]]
      args  <- as.list(rows[i, params, drop = FALSE])
      args  <- args[!vapply(args, is.na, logical(1))]

      component <- function(h) do.call(model, c(list(h), args))
      component(numeric(0))          # Parámetros inválidos fallan aquí y no al evaluar
      component
    })

    function(h) Reduce(`+`, lapply(structures, function(component) component(h)))
  })

  stats::setNames(variograms, ids)
}
