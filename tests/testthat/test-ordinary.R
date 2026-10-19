sim_data <- function(sites, seed, n = 50) {
  set.seed(seed)
  t <- seq_len(n)
  dat <- sapply(sites, function(i) {
    rnorm(1, 5) * sin(2 * pi * t / 50) + rnorm(1) * cos(2 * pi * t / 50) + rnorm(n, sd = 0.1)
  })
  colnames(dat) <- sites
  dat
}

sites  <- paste0("s", 1:10)
set.seed(1)
coords <- cbind(x = runif(10), y = runif(10))
rownames(coords) <- sites

a  <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")
b  <- spatfd(sim_data(sites[1:8], 2), coords, basis = "fourier", nbasis = 11, nharm = 1, name = "b")
ab <- a + b

# Tres sitios nuevos y uno observado (s3)
new <- rbind(cbind(x = c(0.2, 0.7, 0.5), y = c(0.5, 0.5, 0.1)), coords["s3", , drop = FALSE])
rownames(new) <- c("n1", "n2", "n3", "s3")

fitted <- list(
  a.PC1 = gstat::vgm(30, "Exp", 0.5, nugget = 2),
  a.PC2 = gstat::vgm(1, "Sph", 0.8)
)
a <- add_model(a, fitted, name = "single")

# Datos en el formato de gstat, sin los sitios con score faltante
gstat_data <- function(sfd, id) {
  d <- as.data.frame(sfd$scores)
  d <- d[!is.na(d[[id]]), ]
  sp::coordinates(d) <- ~ x + y
  d
}
gstat_new <- as.data.frame(new)
sp::coordinates(gstat_new) <- ~ x + y

# Modelo lineal de corregionalización de gstat. Sin 'beta', gstat hace kriging
# (y cokriging) ordinario.
gstat_lmc <- function(sfd) {
  g <- NULL
  for (id in c("a.PC1", "a.PC2", "b.PC1")) {
    g <- gstat::gstat(g, id, reformulate("1", id), gstat_data(sfd, id))
  }
  g <- gstat::gstat(g, id = "a.PC1", model = gstat::vgm(30, "Exp", 0.5, nugget = 2))
  g <- gstat::gstat(g, id = "a.PC2", model = gstat::vgm(3, "Exp", 0.5, nugget = 1))
  g <- gstat::gstat(g, id = "b.PC1", model = gstat::vgm(20, "Exp", 0.5, nugget = 2))
  g <- gstat::gstat(g, id = c("a.PC1", "a.PC2"), model = gstat::vgm(2, "Exp", 0.5, nugget = 0.5))
  g <- gstat::gstat(g, id = c("a.PC1", "b.PC1"), model = gstat::vgm(-10, "Exp", 0.5, nugget = 1))
  gstat::gstat(g, id = c("a.PC2", "b.PC1"), model = gstat::vgm(1, "Exp", 0.5, nugget = -0.5))
}

# Curvas observadas evaluadas en los tiempos t, sin centrar y su media
observed_curves <- function(sfd, v, t = 1:50) {
  observed <- fda::eval.fd(t, get_fd(sfd, v))
  colnames(observed) <- colnames(get_fd(sfd, v)$coefs)
  observed
}
mean_curve <- function(sfd, v, t = 1:50) as.vector(fda::eval.fd(t, get_pca(sfd, v)$meanfd))

gstat_covariance <- function(m, d) gstat::variogramLine(m, dist_vector = d, covariance = TRUE)
distances <- function(p, q) {
  as.matrix(dist(rbind(p, q)))[seq_len(nrow(p)), nrow(p) + seq_len(nrow(q)), drop = FALSE]
}

# ................................ MÉTODO SCORES ................................

test_that("ordinary kriging is the default and reproduces the ordinary kriging of gstat", {
  kriging <- krige_scores(a, new, model = "single")

  for (id in names(fitted)) {
    ref <- gstat::krige(reformulate("1", id), gstat_data(a, id), gstat_new,
                        model = fitted[[id]], debug.level = 0)
    expect_equal(unname(kriging$pred[, id]), ref$var1.pred, info = id)
    expect_equal(unname(kriging$variance[, id]), ref$var1.var, info = id)
  }

  pred <- predict(a, new, model = "single")
  expect_equal(attr(pred, "kriging"), "ordinary")
  expect_equal(pred$a$var, rowSums(kriging$variance))
  expect_output(print(pred), "scores, ordinary kriging")

  # Es distinto del kriging simple, salvo en el sitio observado
  simple <- predict(a, new, model = "single", kriging = "simple")
  expect_false(isTRUE(all.equal(pred$a$pred[, "n1"], simple$a$pred[, "n1"])))
  expect_equal(pred$a$pred[, "s3"], simple$a$pred[, "s3"])
  expect_equal(pred$a$var[["s3"]], 0)
})

test_that("ordinary cokriging reproduces the one of gstat", {
  g      <- gstat_lmc(ab)
  ref    <- as.data.frame(predict(g, gstat_new, debug.level = 0))
  ab_lmc <- add_model(ab, toSpatFD(g), name = "lmc")

  kriging <- krige_scores(ab_lmc, new, model = "lmc")
  for (id in c("a.PC1", "a.PC2", "b.PC1")) {
    expect_equal(unname(kriging$pred[, id]), ref[[paste0(id, ".pred")]], info = id)
    expect_equal(unname(kriging$variance[, id]), ref[[paste0(id, ".var")]], info = id)
  }
  expect_true(attr(predict(ab_lmc, new, model = "lmc"), "crossed"))
})

test_that("models without a sill are used with ordinary kriging only", {
  unbounded <- list(a.PC1 = gstat::vgm(20, "Lin", 0, nugget = 1),
                    a.PC2 = gstat::vgm(2, "Pow", 1.5))
  sfd <- add_model(a, unbounded, name = "unbounded")

  kriging <- krige_scores(sfd, new, model = "unbounded")
  for (id in names(unbounded)) {
    ref <- gstat::krige(reformulate("1", id), gstat_data(sfd, id), gstat_new,
                        model = unbounded[[id]], debug.level = 0)
    expect_equal(unname(kriging$pred[, id]), ref$var1.pred, info = id)
    expect_equal(unname(kriging$variance[, id]), ref$var1.var, info = id)
  }

  expect_error(predict(sfd, new, model = "unbounded", kriging = "simple"),
               "has no sill.*Use kriging = 'ordinary'")
  expect_error(predict(sfd, new, kriging = "other"), "should be one of")
})

# ................................ MÉTODO LAMBDA ................................

test_that("the ordinary lambda method does not depend on centring the curves", {
  pred <- predict(a, new, model = "single", method = "lambda")
  expect_equal(attr(pred, "kriging"), "ordinary")

  # Referencia: el sistema con la restricción de que los pesos sumen 1, escrito
  # directamente con las covarianzas de gstat
  places <- a$scores[, c("x", "y")]
  n      <- nrow(places)
  Omega  <- Reduce(`+`, lapply(fitted, gstat_covariance, d = distances(places, places)))
  C0     <- Reduce(`+`, lapply(fitted, gstat_covariance, d = distances(places, new)))
  sill   <- sum(vapply(fitted, function(m) sum(m$psill), numeric(1)))

  system <- rbind(cbind(Omega, 1), c(rep(1, n), 0))
  sol    <- solve(system, rbind(C0, 1))
  lambda <- sol[seq_len(n), ]
  mu     <- sol[n + 1, ]

  expect_equal(unname(colSums(lambda)), rep(1, nrow(new)))
  expect_equal(unname(pred$a$var), unname(sill - colSums(C0 * lambda) - mu))

  # Centradas o no, la curva es la misma: la media se cancela porque los pesos
  # suman 1
  raw     <- observed_curves(a, "a")[, rownames(places)]
  centred <- mean_curve(a, "a") + (raw - mean_curve(a, "a")) %*% lambda
  expect_equal(unname(pred$a$pred), unname(raw %*% lambda))
  expect_equal(unname(pred$a$pred), unname(centred))

  expect_true(all(pred$a$var >= -1e-8))
  expect_equal(pred$a$var[["s3"]], 0, tolerance = 1e-8)
  expect_equal(unname(pred$a$pred[, "s3"]), unname(raw[, "s3"]))
})

test_that("the ordinary multivariate lambda method solves the constrained cokriging system", {
  lmc    <- gstat_lmc(ab)
  ab_lmc <- add_model(ab, toSpatFD(lmc), name = "lmc")
  pred   <- predict(ab_lmc, new, model = "lmc", method = "lambda")

  G <- inner_product(get_pca(ab, "a")$harmonics, get_pca(ab, "b")$harmonics, c(1, 50))

  sites_a <- ab$scores[!is.na(ab$scores[, "a.PC1"]), c("x", "y")]
  sites_b <- ab$scores[!is.na(ab$scores[, "b.PC1"]), c("x", "y")]
  C <- function(id, d) gstat_covariance(lmc$model[[id]], d)
  cov_aa <- function(d) C("a.PC1", d) + C("a.PC2", d)
  cov_bb <- function(d) C("b.PC1", d)
  cov_ab <- function(d) G[1, 1] * C("a.PC1.b.PC1", d) + G[2, 1] * C("a.PC2.b.PC1", d)

  na <- nrow(sites_a)
  nb <- nrow(sites_b)
  Sigma <- rbind(cbind(cov_aa(distances(sites_a, sites_a)),    cov_ab(distances(sites_a, sites_b))),
                 cbind(t(cov_ab(distances(sites_a, sites_b))), cov_bb(distances(sites_b, sites_b))))
  F <- rbind(cbind(rep(1, na), 0), cbind(0, rep(1, nb)))
  system <- rbind(cbind(Sigma, F), cbind(t(F), matrix(0, 2, 2)))

  raw <- cbind(observed_curves(ab, "a")[, rownames(sites_a)], observed_curves(ab, "b")[, rownames(sites_b)])
  m   <- nrow(new)
  targets <- list(
    a = list(C0 = rbind(cov_aa(distances(sites_a, new)), cov_ab(distances(sites_b, new))),
             e = rbind(rep(1, m), rep(0, m)), sill = cov_aa(matrix(0))),
    b = list(C0 = rbind(cov_ab(distances(sites_a, new)), cov_bb(distances(sites_b, new))),
             e = rbind(rep(0, m), rep(1, m)), sill = cov_bb(matrix(0)))
  )
  for (v in names(targets)) {
    rhs    <- rbind(targets[[v]]$C0, targets[[v]]$e)
    sol    <- solve(system, rhs)
    lambda <- sol[seq_len(na + nb), ]

    # Los pesos de la variable que se predice suman 1 y los de la otra 0, así
    # que la predicción usa las curvas sin centrar
    own <- if (v == "a") seq_len(na) else na + seq_len(nb)
    expect_equal(unname(colSums(lambda[own, ])), rep(1, m), info = v)
    expect_equal(unname(colSums(lambda[-own, ])), rep(0, m), info = v)

    expect_equal(unname(pred[[v]]$pred), unname(raw %*% lambda), info = v)
    expect_equal(unname(pred[[v]]$var), unname(drop(targets[[v]]$sill) - colSums(rhs * sol)), info = v)
  }
})

test_that("the lambda method also works with models without a sill", {
  sfd  <- add_model(a, list(a.PC1 = gstat::vgm(20, "Lin", 0, nugget = 1),
                            a.PC2 = gstat::vgm(2, "Pow", 1.5)), name = "unbounded")
  pred <- predict(sfd, new, model = "unbounded", method = "lambda")

  expect_true(all(pred$a$var >= -1e-8))
  expect_equal(pred$a$var[["s3"]], 0, tolerance = 1e-8)
  expect_equal(unname(pred$a$pred[, "s3"]), unname(observed_curves(a, "a")[, "s3"]))
  expect_error(predict(sfd, new, model = "unbounded", method = "lambda", kriging = "simple"),
               "has no sill")
})

# ............................ VALIDACIÓN CRUZADA ..............................

test_that("leave_one_out() uses ordinary kriging by default", {
  sfd <- add_model(a, list(a.PC1 = gstat::vgm(20, "Lin", 0, nugget = 1),
                           a.PC2 = gstat::vgm(2, "Pow", 1.5)), name = "unbounded")

  loo <- leave_one_out(sfd)
  expect_equal(attr(loo, "kriging"), "ordinary")
  expect_setequal(paste(loo$summary$model, loo$summary$method),
                  c("single scores", "single lambda", "unbounded scores", "unbounded lambda"))
  expect_output(print(loo), "ordinary kriging")

  # Con kriging simple quedan por fuera los modelos sin meseta
  simple <- leave_one_out(sfd, kriging = "simple")
  expect_setequal(simple$summary$model, "single")
  expect_false(isTRUE(all.equal(simple$errors[, "single / scores"], loo$errors[, "single / scores"])))
})
