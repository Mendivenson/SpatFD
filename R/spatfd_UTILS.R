# ======================= FUNCIONES INTERNAS DE spatfd() =======================
# Solo las usa spatfd(). Lo compartido por varias funciones está en utils.R.

# SUAVIZAMIENTO CON DATOS FALTANTES: fda::smooth.basis() no falla con NA, pero
# devuelve todos los coeficientes de esa curva como NA. Las curvas completas se
# suavizan juntas y las que tienen faltantes una por una, solo con sus tiempos
# observados (mínimos cuadrados penalizados no necesita que todas las curvas
# compartan tiempos).
smooth_curves <- function(t, dat, fdpar){
  complete <- colSums(is.na(dat)) == 0
  coefs    <- matrix(NA_real_, fdpar$fd$basis$nbasis, ncol(dat), dimnames = list(NULL, colnames(dat)))

  if (any(complete)) {
    coefs[, complete] <- fda::smooth.basis(t, dat[, complete, drop = F], fdpar)$fd$coefs
  }
  for (j in which(!complete)) {
    observed   <- !is.na(dat[, j])
    coefs[, j] <- fda::smooth.basis(t[observed], dat[observed, j], fdpar)$fd$coefs
  }

  fda::fd(coefs, fdpar$fd$basis)
}
