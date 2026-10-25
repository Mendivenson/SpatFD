# ====================== FUNCIONES INTERNAS DE to_gstat() ======================
# Solo las usa to_gstat(). Lo compartido por varias funciones está en utils.R.

# NÚMERO DE PC POR VARIABLE: Vector con nombres, uno por variable, completado
# con todas las PC de las variables que no aparecen en 'pcs'
check_pcs <- function(objects, pcs){
  available <- vapply(objects, function(sfd) ncol(get_pca(sfd)$scores), integer(1))
  if (is.null(pcs)) return(available)

  stopifnot(
    "'pcs' must be a named numeric vector" = is.numeric(pcs) && !is.null(names(pcs)) && all(names(pcs) != '')
  )

  unknown <- setdiff(names(pcs), names(objects))
  if (length(unknown) > 0) {
    stop("Variables in 'pcs' not given in '...': ", paste(unknown, collapse = ', '), call. = F)
  }

  invalid <- names(pcs)[pcs < 1 | pcs != round(pcs) | pcs > available[names(pcs)]]
  if (length(invalid) > 0) {
    stop("'pcs' must be a whole number between 1 and the number of components of the variable. Invalid for: ",
         paste0(invalid, ' (', available[invalid], ' components)', collapse = ', '), call. = F)
  }

  available[names(pcs)] <- as.integer(pcs)
  available
}
