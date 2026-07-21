togstat <- function(sfd, vars = var_names(sfd)) {
  stopifnot(
    'function only available for SpatFD objects' = inherits(sfd, 'SpatFD'),
    'Variables not found in object' = all(vars %in% var_names(sfd))
  )

  # SEPARAR COORDS Y SCORES:
  scores <- sfd$scores[, !(colnames(sfd$scores) %in% colnames(attr(sfd, 'data')$coords))]
  coords <- sfd$scores[, colnames(attr(sfd, 'data')$coords)]

  # Selección de los componentes principales con los que se va a trabajar
  scores <- scores[, grepl(paste0("^(", paste0(vars, '\\.', collapse = '|'), ")"),
                           colnames(scores)), drop = FALSE]

  # Añadir coordenadas a los scores
  scores <- as.data.frame(scores)
  sp::coordinates(scores) <- coords[rownames(scores), ]

  # Generar el objeto gstat
  g <- NULL
  for (j in names(scores)) {
    g <- gstat::gstat(g       = g,
                      id      = j,
                      formula = reformulate('1', response = j),
                      data    = scores[!is.na(scores[[j]]), j])
  }

  # Generar una versión básica del variograma
  # vg <- gstat::variogram(g, cross = length(vars) != 1)

  return(g)
}
