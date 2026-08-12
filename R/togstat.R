togstat <- function(sfd, vars = var_names(sfd)) {

  # CLASS HANDLING
  stopifnot(
    'function only available for SpatFD objects' = inherits(sfd, 'SpatFD'),
    'Variables not found in object' = all(vars %in% var_names(sfd))
  )

  # SEPARAR SCORES Y COORDENADAS:
  # Se separan de esta forma para permitir cualquier nombre de coordenada y variable posible
  scores <- sfd$scores[, !(colnames(sfd$scores) %in% colnames(attr(sfd, 'data')$coords))]
  coords <- attr(sfd, 'data')$coords

  # SELECCIÓN DE LOS COMPONENTES PRINCIPALES CON LOS QUE SE VA A TRABAJAR:
  # Básicamente se seleccionan las variables escogidas por el usuario.
  scores <- scores[, grepl(paste0("^(", paste0(vars, '\\.', collapse = '|'), ")"),
                           colnames(scores)), drop = FALSE]

  # SE AÑADEN LAS COORDENADAS ESPACIALES A LOS SCORES.
  scores <- as.data.frame(scores)
  sp::coordinates(scores) <- coords[rownames(scores), ]

  # GENERACIÓN OBJETO GSTAT:
  # El objeto gstat se genera en conjunto independientemente de si se seleccionan una o varias variables
  g <- NULL
  for (j in names(scores)) {
    g <- gstat::gstat(g       = g,
                      id      = j,
                      formula = reformulate('1', response = j),
                      data    = scores[!is.na(scores[[j]]), j])
  }

  return(g)
}
