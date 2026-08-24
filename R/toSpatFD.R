toSpatFD <- function(model, ...) {
  UseMethod("toSpatFD")
}

# MÉTODO IMPLEMENTADO PARA LOS MODELOS PROVENIENTES DE GSTAT.
toSpatFD.gstat <- function(model, sfd, ...) {
  model.trans <- c(
    'Nug' = 'nugget',
    'Exp' = 'exponential',
    'Sph' = 'spherical',
    'Gau' = 'gaussian',
    'Exc' = 'stable',
    'Mat' = 'matern',
    'Ste' = "stein's matern",
    'Cir' = 'circular',
    'Lin' = 'linear',
    'Bes' = 'bessel',
    'Pen' = 'pentamodel',
    'Per' = 'periodic',
    'Wav' = 'wave',
    'Hol' = 'hole',
    'Log' = 'logarithmic',
    'Pow' = 'power',
    'Spl' = 'spline',
    'Leg' = 'legendre'
  )

  stopifnot(
    'There are no models to extract from this object' = length(model$model) != 0
  )

  model <- model$model
  vars  <- rep(names(model), each = unlist(lapply(model, nrow)))
  model <- model |>
    do.call(what = rbind)
  model <- cbind('var' = rownames(model), model)
  model <- model[,c('var', 'model', 'range', 'psill', 'kappa')]
  colnames(model) <- c('var', 'model', 'rango', 'sigma', 'alpha1')
  rownames(model) <- NULL

  model[,"model"] = model.trans[model$model]
  model$var = vars

  return(model)
}

# CUANDO SE TENGAN LISTAS DE MODELOS EN EL CASO UNIVARIADO
toSpatFD.list <- function(model,...){

  # Modelos GSTAT
  if (all(c('variogramModel', 'data.frame') == lapply(model, class) |>
      unlist() |>
      unique())){
        x = list()
        x$model = model
        class(x) = c('gstat', class(x))
  }

  toSpatFD(x)
}

toSpatFD.default <- function(model, ...) {
  stop("No SpatFD translator registered for objects of class '",
       class(model)[1], "'", call. = FALSE)
}
