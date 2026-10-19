sim_data <- function(sites, seed) {
  set.seed(seed)
  t <- 1:50
  dat <- sapply(sites, function(i) {
    rnorm(1, 5) * sin(2 * pi * t / 50) + rnorm(1) * cos(2 * pi * t / 50) + rnorm(50, sd = 0.1)
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

# Datos en el formato de gstat, sin los sitios con score faltante
gstat_data <- function(sfd, id) {
  d <- as.data.frame(sfd$scores)
  d <- d[!is.na(d[[id]]), ]
  sp::coordinates(d) <- ~ x + y
  d
}
gstat_new <- as.data.frame(new)
sp::coordinates(gstat_new) <- ~ x + y

a <- add_model(a, toSpatFD(fitted), name = "single")

test_that("univariate prediction reproduces the simple kriging of gstat", {
  kriging <- krige_scores(a, new, model = "single", kriging = "simple")

  for (id in names(fitted)) {
    ref <- gstat::krige(reformulate("1", id), gstat_data(a, id), gstat_new,
                        model = fitted[[id]], beta = 0, debug.level = 0)
    expect_equal(unname(kriging$pred[, id]), ref$var1.pred, info = id)
    expect_equal(unname(kriging$variance[, id]), ref$var1.var, info = id)
  }

  # La varianza que ve el usuario es la suma de las varianzas de los scores
  expect_equal(predict(a, new, model = "single", kriging = "simple")$a$var, rowSums(kriging$variance))
})

# Modelo lineal de corregionalización de gstat para los tres scores de ab, con
# media cero (kriging simple)
gstat_lmc <- function(sfd) {
  g <- NULL
  for (id in c("a.PC1", "a.PC2", "b.PC1")) {
    g <- gstat::gstat(g, id, reformulate("1", id), gstat_data(sfd, id), beta = 0)
  }
  g <- gstat::gstat(g, id = "a.PC1", model = gstat::vgm(30, "Exp", 0.5, nugget = 2))
  g <- gstat::gstat(g, id = "a.PC2", model = gstat::vgm(3, "Exp", 0.5, nugget = 1))
  g <- gstat::gstat(g, id = "b.PC1", model = gstat::vgm(20, "Exp", 0.5, nugget = 2))
  g <- gstat::gstat(g, id = c("a.PC1", "a.PC2"), model = gstat::vgm(2, "Exp", 0.5, nugget = 0.5))
  g <- gstat::gstat(g, id = c("a.PC1", "b.PC1"), model = gstat::vgm(-10, "Exp", 0.5, nugget = 1))
  gstat::gstat(g, id = c("a.PC2", "b.PC1"), model = gstat::vgm(1, "Exp", 0.5, nugget = -0.5))
}

test_that("multivariate prediction reproduces the simple cokriging of gstat", {
  g    <- gstat_lmc(ab)
  ref  <- as.data.frame(predict(g, gstat_new, debug.level = 0))
  ab_lmc  <- add_model(ab, toSpatFD(g), name = "lmc")
  kriging <- krige_scores(ab_lmc, new, model = "lmc", kriging = "simple")

  for (id in c("a.PC1", "a.PC2", "b.PC1")) {
    expect_equal(unname(kriging$pred[, id]), ref[[paste0(id, ".pred")]], info = id)
    expect_equal(unname(kriging$variance[, id]), ref[[paste0(id, ".var")]], info = id)
  }

  pred <- predict(ab_lmc, new, model = "lmc", kriging = "simple")
  expect_equal(pred$a$var, rowSums(kriging$variance[, c("a.PC1", "a.PC2")]))
  expect_equal(pred$b$var, kriging$variance[, "b.PC1"])
  expect_true(attr(pred, "crossed"))
})

test_that("without cross-variograms the scores are predicted separately", {
  ab_direct <- add_model(ab, toSpatFD(c(fitted, list(b.PC1 = gstat::vgm(20, "Exp", 0.5)))),
                         name = "direct")

  joint <- predict(ab_direct, new, model = "direct")
  alone <- predict(ab_direct, new, model = "direct", variable = "a")

  expect_named(joint, c("a", "b"))
  expect_named(alone, "a")

  expect_equal(joint$a$pred, alone$a$pred)
  expect_equal(joint$a$var, alone$a$var)
})

test_that("observed sites are predicted exactly", {
  kriging <- krige_scores(a, new, model = "single")
  pred    <- predict(a, new, model = "single")$a

  expect_equal(kriging$pred["s3", ], a$scores["s3", c("a.PC1", "a.PC2")])
  expect_equal(pred$var[["s3"]], 0)
  expect_equal(unname(pred$pred[, "s3"]),
               unname(fda::eval.fd(1:50, get_fd(a, "a"))[, "s3"]),
               tolerance = 0.1)
})

test_that("the predicted curves are matrices rebuilt from the predicted scores", {
  pred <- predict(a, new, model = "single")
  pca  <- get_pca(a, "a")
  t    <- 1:50

  expected <- fda::eval.fd(t, pca$harmonics) %*% t(krige_scores(a, new, model = "single")$pred)
  expected <- expected + as.vector(fda::eval.fd(t, pca$meanfd))

  # Queda registrado cómo se predijo
  expect_equal(attr(pred, "method"), "scores")
  expect_false(attr(pred, "crossed"))

  # Siempre una lista por variable con las curvas (pred) y su varianza (var).
  # Las curvas son una matriz tiempos x sitios, como los datos de entrada.
  expect_s3_class(pred, "SpatFDpred")
  expect_named(pred, "a")
  expect_named(pred$a, c("pred", "var"))
  expect_true(is.matrix(pred$a$pred))
  expect_equal(dim(pred$a$pred), c(length(t), nrow(new)))
  expect_equal(dimnames(pred$a$pred), list(as.character(t), rownames(new)))
  expect_named(pred$a$var, rownames(new))

  expect_equal(unname(pred$a$pred), unname(expected))

  # Las curvas funcionales solo quedan en el atributo oculto
  expect_setequal(names(attributes(pred)),
                  c("names", "class", "newcoords", "method", "kriging", "crossed", "internal"))
  expect_s3_class(predicted_fd(pred, "a"), "fd")
})

test_that("predict() validates its arguments", {
  sfd <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")

  expect_error(predict(sfd, new), "has no models")

  # Con un único modelo registrado no hace falta nombrarlo
  sfd <- add_model(sfd, toSpatFD(fitted), name = "single")
  expect_equal(predict(sfd, new)$a$pred, predict(sfd, new, model = "single")$a$pred)

  expect_error(predict(sfd, new, variable = "z"), "Variables not found")
  expect_error(predict(sfd, new, method = c("scores", "lambda")), "must be a single string")
  expect_error(predict(sfd, new[, "x", drop = FALSE]), "must have the coordinates of the object")

  # En los errores sobre el modelo se muestra la tabla resumen
  sfd <- add_model(sfd, toSpatFD(fitted["a.PC1"]), name = "incomplete")
  expect_error(predict(sfd, new), "several models.*single.*incomplete")
  expect_error(predict(sfd, new, model = "other"), "must be the name of one of the models.*single")
  expect_error(predict(sfd, new, model = "incomplete", variable = "a"),
               "no variogram for the scores: a.PC2.*single")

  # El kriging simple necesita modelos con meseta
  sfd <- add_model(sfd, toSpatFD(list(a.PC1 = gstat::vgm(2, "Lin", 0), a.PC2 = fitted$a.PC2)),
                   name = "unbounded")
  expect_error(predict(sfd, new, model = "unbounded", kriging = "simple"), "has no sill")
})

test_that("add_model() returns the object with the model registered", {
  sfd <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")
  expect_equal(nrow(sfd$models), 0)

  # Como cualquier función de modelado de R: no modifica el objeto original
  with_model <- add_model(sfd, toSpatFD(fitted))
  expect_equal(nrow(sfd$models), 0)
  expect_equal(nrow(with_model$models), 1)

  # Nombres secuenciales por defecto, o el que dé el usuario
  sfd <- sfd |>
    add_model(toSpatFD(fitted)) |>
    add_model(toSpatFD(fitted["a.PC1"]), name = "only PC1") |>
    add_model(toSpatFD(fitted))

  expect_equal(rownames(sfd$models), c("model 1", "only PC1", "model 3"))
  expect_named(sfd$models, c("variables", "structures", "no. scores", "cross", "sill"))
  expect_equal(sfd$models$structures, c("nugget, exponential, spherical", "nugget, exponential",
                                        "nugget, exponential, spherical"))
  expect_equal(sfd$models[["no. scores"]], c(2, 1, 2))
  expect_equal(sfd$models$sill, c(33, 32, 33))
  expect_equal(sfd[["models"]], sfd$models)

  # Los modelos compilados quedan en el atributo oculto
  compiled <- get_models(sfd)[["model 1"]]
  expect_named(compiled$variograms, c("a.PC1", "a.PC2"))
  expect_equal(compiled$variograms$a.PC1(c(0, 1e6)), c(0, 32))

  # Un nombre repetido reemplaza el modelo sin cambiar su posición
  expect_message(sfd <- add_model(sfd, toSpatFD(fitted["a.PC2"]), name = "only PC1"), "Replacing")
  expect_equal(rownames(sfd$models), c("model 1", "only PC1", "model 3"))
  expect_equal(sfd$models$structures[2], "spherical")
})

test_that("add_model() translates fitted gstat models on the way in", {
  sfd <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")

  sfd <- add_model(sfd, fitted, name = "gstat list")
  sfd <- add_model(sfd, toSpatFD(fitted), name = "table")

  expect_equal(get_models(sfd)[["gstat list"]]$table, get_models(sfd)[["table"]]$table)
  expect_equal(predict(sfd, new, model = "gstat list")$a$pred,
               predict(sfd, new, model = "table")$a$pred)

  expect_error(add_model(sfd, 1), "No SpatFD translator")
})

test_that("add_model() validates the models", {
  sfd <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")

  expect_error(add_model(toSpatFD(fitted)), "only available for SpatFD objects")
  expect_error(add_model(sfd, toSpatFD(list(z.PC1 = fitted$a.PC1))), "do not match the scores")

  df <- toSpatFD(fitted)
  df$model[1] <- "other"
  expect_error(add_model(sfd, df), "Unknown variogram models: other")
})

test_that("models are kept when objects are combined", {
  x <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")
  y <- spatfd(sim_data(sites[1:8], 2), coords, basis = "fourier", nbasis = 11, nharm = 1, name = "b")
  x <- add_model(x, toSpatFD(fitted))
  y <- add_model(y, toSpatFD(list(b.PC1 = gstat::vgm(20, "Exp", 0.5))))

  expect_message(xy <- x + y, "renamed to 'model 2'")
  expect_equal(rownames(xy$models), c("model 1", "model 2"))
  expect_equal(xy$models$variables, c("a", "b"))
  expect_false(any(xy$models$cross))

  # El objeto combinado es independiente de los originales
  xy <- add_model(xy, toSpatFD(fitted))
  expect_equal(nrow(xy$models), 3)
  expect_equal(nrow(x$models), 1)
})

test_that("the nugget variogram is zero at the origin", {
  expect_equal(models$nugget(c(0, 1, 5), sigma = 2), c(0, 2, 2))
})

test_that("predictions are plotted over the observed curves", {
  pred <- predict(a, new, model = "single")

  grDevices::pdf(NULL)
  expect_invisible(plot(pred))
  expect_no_error(plot(pred, variable = "a", col = "red", lwd = 2))
  grDevices::dev.off()

  expect_error(plot(pred, variable = "z"), "was not predicted")
  expect_error(plot(pred, mode = "other"), "'mode' must be")
  expect_error(plot(pred, mode = "map"), "only available for predictions made over a map")
})

test_that("a map is predicted on a regular grid and drawn as two maps", {
  skip_if_not_installed("sf")
  square <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))))

  # predict() no dibuja: devuelve la predicción en cerca de n puntos dentro del mapa
  pred <- predict(a, square, model = "single", n = 100)
  grid <- attr(pred, "newcoords")

  expect_s3_class(pred, "SpatFDpred")
  expect_equal(nrow(grid), 100, tolerance = 0.2)
  expect_equal(colnames(grid), c("x", "y"))
  expect_true(all(grid >= 0 & grid <= 1))
  expect_length(pred$a$var, nrow(grid))
  expect_equal(ncol(pred$a$pred), nrow(grid))

  # Es la misma predicción que se obtiene dando esos puntos directamente
  direct <- predict(a, grid, model = "single")
  expect_equal(pred$a$pred, direct$a$pred)
  expect_equal(pred$a$var, direct$a$var)

  grDevices::pdf(NULL)
  expect_invisible(plot(pred))                                       # mapas por defecto
  expect_no_error(plot(pred, time = 25))
  expect_no_error(plot(pred, time = 25, main = "My title"))
  expect_no_error(plot(pred, mode = "curves"))
  expect_no_error(plot(pred, mode = "map", col = grDevices::heat.colors(20)))
  grDevices::dev.off()

  expect_error(plot(pred, time = 1000), "'time' must be a single value within the range")
  expect_error(predict(a, square, model = "single", n = -1), "'n' must be a single positive number")
})

# ................................ MÉTODO LAMBDA ................................

# Curvas observadas de una variable, centradas y evaluadas en los tiempos t
centred_curves <- function(sfd, v, t = 1:50) {
  observed <- fda::eval.fd(t, get_fd(sfd, v))
  colnames(observed) <- colnames(get_fd(sfd, v)$coefs)
  observed - as.vector(fda::eval.fd(t, get_pca(sfd, v)$meanfd))
}

# Covarianza de un modelo de gstat en una matriz de distancias
gstat_covariance <- function(m, d) gstat::variogramLine(m, dist_vector = d, covariance = TRUE)

test_that("the univariate lambda method is simple kriging of the centred curves", {
  pred <- predict(a, new, model = "single", method = "lambda", kriging = "simple")

  expect_equal(attr(pred, "method"), "lambda")
  expect_false(attr(pred, "crossed"))
  expect_named(pred, "a")
  expect_named(pred$a, c("pred", "var"))
  expect_equal(dimnames(pred$a$pred), list(as.character(1:50), rownames(new)))

  # Referencia: las fórmulas escritas directamente, con las covarianzas de gstat
  places <- a$scores[, c("x", "y")]
  d      <- as.matrix(dist(rbind(places, new)))
  d_obs  <- d[seq_len(nrow(places)), seq_len(nrow(places))]
  d_new  <- d[seq_len(nrow(places)), -seq_len(nrow(places)), drop = FALSE]

  Omega  <- Reduce(`+`, lapply(fitted, gstat_covariance, d = d_obs))
  C0     <- Reduce(`+`, lapply(fitted, gstat_covariance, d = d_new))
  lambda <- solve(Omega, C0)

  expected <- as.vector(fda::eval.fd(1:50, get_pca(a, "a")$meanfd)) +
    centred_curves(a, "a")[, rownames(places)] %*% lambda
  sill <- sum(vapply(fitted, function(m) sum(m$psill), numeric(1)))

  expect_equal(unname(pred$a$pred), unname(expected))
  expect_equal(unname(pred$a$var), unname(sill - colSums(C0 * lambda)))

  # Con la meseta del modelo la varianza no puede ser negativa, y en un sitio
  # observado la curva predicha es la observada
  expect_true(all(pred$a$var >= -1e-8))
  expect_equal(pred$a$var[["s3"]], 0, tolerance = 1e-8)
  expect_equal(unname(pred$a$pred[, "s3"]), unname(fda::eval.fd(1:50, get_fd(a, "a"))[, 3]))
})

test_that("the multivariate lambda method solves the functional cokriging system", {
  lmc    <- gstat_lmc(ab)
  ab_lmc <- add_model(ab, toSpatFD(lmc), name = "lmc")
  pred   <- predict(ab_lmc, new, model = "lmc", method = "lambda", kriging = "simple")

  expect_named(pred, c("a", "b"))
  expect_true(attr(pred, "crossed"))

  # Referencia: el sistema por bloques de Bohorquez, Giraldo y Mateu (2017),
  # con G el producto interno entre las funciones propias de las dos
  # variables. Dentro de una variable G es la identidad, así que el cruzado
  # a.PC1.a.PC2 no interviene.
  G <- inner_product(get_pca(ab, "a")$harmonics, get_pca(ab, "b")$harmonics, c(1, 50))
  expect_equal(inner_product(get_pca(ab, "a")$harmonics, get_pca(ab, "a")$harmonics, c(1, 50)),
               diag(2), ignore_attr = TRUE)

  sites_a <- ab$scores[!is.na(ab$scores[, "a.PC1"]), c("x", "y")]
  sites_b <- ab$scores[!is.na(ab$scores[, "b.PC1"]), c("x", "y")]
  distances <- function(p, q) {
    as.matrix(dist(rbind(p, q)))[seq_len(nrow(p)), nrow(p) + seq_len(nrow(q)), drop = FALSE]
  }
  C <- function(id, d) gstat_covariance(lmc$model[[id]], d)

  cov_aa <- function(d) C("a.PC1", d) + C("a.PC2", d)
  cov_bb <- function(d) C("b.PC1", d)
  cov_ab <- function(d) G[1, 1] * C("a.PC1.b.PC1", d) + G[2, 1] * C("a.PC2.b.PC1", d)

  Sigma <- rbind(cbind(cov_aa(distances(sites_a, sites_a)),    cov_ab(distances(sites_a, sites_b))),
                 cbind(t(cov_ab(distances(sites_a, sites_b))), cov_bb(distances(sites_b, sites_b))))
  observed <- cbind(centred_curves(ab, "a")[, rownames(sites_a)],
                    centred_curves(ab, "b")[, rownames(sites_b)])

  targets <- list(
    a = list(C0 = rbind(cov_aa(distances(sites_a, new)), cov_ab(distances(sites_b, new))),
             sill = cov_aa(matrix(0))),
    b = list(C0 = rbind(cov_ab(distances(sites_a, new)), cov_bb(distances(sites_b, new))),
             sill = cov_bb(matrix(0)))
  )
  for (v in names(targets)) {
    lambda   <- solve(Sigma, targets[[v]]$C0)
    expected <- as.vector(fda::eval.fd(1:50, get_pca(ab, v)$meanfd)) + observed %*% lambda

    expect_equal(unname(pred[[v]]$pred), unname(expected), info = v)
    expect_equal(unname(pred[[v]]$var), drop(targets[[v]]$sill) - colSums(targets[[v]]$C0 * lambda),
                 info = v)
  }

  # Las dos variables se observaron en s3: ambas curvas se recuperan exactamente
  expect_equal(unname(pred$a$var["s3"]), 0, tolerance = 1e-8)
  expect_equal(unname(pred$b$pred[, "s3"]), unname(fda::eval.fd(1:50, get_fd(ab, "b"))[, 3]))
})

test_that("without cross-variograms between variables lambda predicts each one alone", {
  ab_direct <- add_model(ab, c(fitted, list(b.PC1 = gstat::vgm(20, "Exp", 0.5))), name = "direct")

  joint <- predict(ab_direct, new, model = "direct", method = "lambda")
  alone <- predict(ab_direct, new, model = "direct", method = "lambda", variable = "a")

  expect_false(attr(joint, "crossed"))
  expect_equal(joint$a, alone$a)
  expect_equal(alone$a, predict(a, new, model = "single", method = "lambda")$a)
})

test_that("variables over different time intervals are combined over the common one", {
  short <- spatfd(sim_data(sites[1:8], 2)[1:40, ], coords, t = 1:40,
                  basis = "fourier", nbasis = 11, nharm = 1, name = "b")
  crossed <- toSpatFD(c(fitted, list(b.PC1 = gstat::vgm(20, "Exp", 0.5))))
  crossed <- rbind(crossed, transform(crossed[crossed$var == "b.PC1", ], var = "a.PC1.b.PC1", sigma = 5))
  ab_short <- add_model(a + short, crossed, name = "crossed")

  expect_message(pred <- predict(ab_short, new, model = "crossed", method = "lambda"),
                 "common time interval \\[1, 40\\]")
  expect_equal(rownames(pred$a$pred), as.character(1:40))
  expect_equal(rownames(pred$b$pred), as.character(1:40))
  expect_true(attr(pred, "crossed"))
  expect_equal(pred$a$var[["s3"]], 0, tolerance = 1e-8)

  # El método scores no se ve afectado
  expect_equal(rownames(predict(ab_short, new, model = "crossed")$a$pred), as.character(1:50))

  grDevices::pdf(NULL)
  expect_no_error(plot(pred))
  grDevices::dev.off()

  # Sin tiempos en común no se pueden combinar
  later <- spatfd(sim_data(sites[1:8], 2), coords, t = 101:150,
                  basis = "fourier", nbasis = 11, nharm = 1, name = "b")
  ab_later <- add_model(a + later, crossed, name = "crossed")
  expect_error(predict(ab_later, new, model = "crossed", method = "lambda"), "no time in common")
})
