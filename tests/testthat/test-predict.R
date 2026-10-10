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

# Los modelos se registran en el sitio, sin reasignar el objeto
add_model(a, toSpatFD(fitted), name = "single")

test_that("univariate prediction reproduces the simple kriging of gstat", {
  kriging <- krige_scores(a, new, model = "single")

  for (id in names(fitted)) {
    ref <- gstat::krige(reformulate("1", id), gstat_data(a, id), gstat_new,
                        model = fitted[[id]], beta = 0, debug.level = 0)
    expect_equal(unname(kriging$pred[, id]), ref$var1.pred, info = id)
    expect_equal(unname(kriging$variance[, id]), ref$var1.var, info = id)
  }

  # La varianza que ve el usuario es la suma de las varianzas de los scores
  expect_equal(predict(a, new, model = "single")$a$var, rowSums(kriging$variance))
})

test_that("multivariate prediction reproduces the simple cokriging of gstat", {
  g <- NULL
  for (id in c("a.PC1", "a.PC2", "b.PC1")) {
    g <- gstat::gstat(g, id, reformulate("1", id), gstat_data(ab, id), beta = 0)
  }
  g <- gstat::gstat(g, id = "a.PC1", model = gstat::vgm(30, "Exp", 0.5, nugget = 2))
  g <- gstat::gstat(g, id = "a.PC2", model = gstat::vgm(3, "Exp", 0.5, nugget = 1))
  g <- gstat::gstat(g, id = "b.PC1", model = gstat::vgm(20, "Exp", 0.5, nugget = 2))
  g <- gstat::gstat(g, id = c("a.PC1", "a.PC2"), model = gstat::vgm(2, "Exp", 0.5, nugget = 0.5))
  g <- gstat::gstat(g, id = c("a.PC1", "b.PC1"), model = gstat::vgm(-10, "Exp", 0.5, nugget = 1))
  g <- gstat::gstat(g, id = c("a.PC2", "b.PC1"), model = gstat::vgm(1, "Exp", 0.5, nugget = -0.5))

  ref  <- as.data.frame(predict(g, gstat_new, debug.level = 0))
  add_model(ab, toSpatFD(g), name = "lmc")
  kriging <- krige_scores(ab, new, model = "lmc")

  for (id in c("a.PC1", "a.PC2", "b.PC1")) {
    expect_equal(unname(kriging$pred[, id]), ref[[paste0(id, ".pred")]], info = id)
    expect_equal(unname(kriging$variance[, id]), ref[[paste0(id, ".var")]], info = id)
  }

  pred <- predict(ab, new, model = "lmc")
  expect_equal(pred$a$var, rowSums(kriging$variance[, c("a.PC1", "a.PC2")]))
  expect_equal(pred$b$var, kriging$variance[, "b.PC1"])
})

test_that("without cross-variograms the scores are predicted separately", {
  add_model(ab, toSpatFD(c(fitted, list(b.PC1 = gstat::vgm(20, "Exp", 0.5)))), name = "direct")

  joint <- predict(ab, new, model = "direct")
  alone <- predict(ab, new, model = "direct", variable = "a")

  expect_named(joint, c("a", "b"))
  expect_named(alone, "a")

  expect_equal(joint$a$pred$coefs, alone$a$pred$coefs)
  expect_equal(joint$a$var, alone$a$var)
})

test_that("observed sites are predicted exactly", {
  kriging <- krige_scores(a, new, model = "single")
  pred    <- predict(a, new, model = "single")$a

  expect_equal(kriging$pred["s3", ], a$scores["s3", c("a.PC1", "a.PC2")])
  expect_equal(pred$var[["s3"]], 0)

  t <- 1:50
  expect_equal(fda::eval.fd(t, pred$pred)[, "s3"],
               fda::eval.fd(t, attr(a, "fd")$a)[, "s3"],
               tolerance = 0.1)
})

test_that("curves are rebuilt from the predicted scores", {
  pred <- predict(a, new, model = "single")
  pca  <- attr(a, "pca")$a
  t    <- 1:50

  expected <- fda::eval.fd(t, pca$harmonics) %*% t(krige_scores(a, new, model = "single")$pred)
  expected <- expected + as.vector(fda::eval.fd(t, pca$meanfd))

  # Siempre una lista por variable con las curvas (pred) y su varianza (var)
  expect_s3_class(pred, "SpatFDpred")
  expect_named(pred, "a")
  expect_named(pred$a, c("pred", "var"))
  expect_s3_class(pred$a$pred, "fd")
  expect_named(pred$a$var, rownames(new))

  expect_equal(unname(fda::eval.fd(t, pred$a$pred)), unname(expected))
  expect_equal(colnames(pred$a$pred$coefs), rownames(new))
})

test_that("predict() validates its arguments", {
  sfd <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")

  expect_error(predict(sfd, new), "has no models")

  # Con un único modelo registrado no hace falta nombrarlo
  add_model(sfd, toSpatFD(fitted), name = "single")
  expect_equal(predict(sfd, new)$a$pred$coefs, predict(sfd, new, model = "single")$a$pred$coefs)

  expect_error(predict(sfd, new, variable = "z"), "Variables not found")
  expect_error(predict(sfd, new, method = "lambda"), "not available yet")
  expect_error(predict(sfd, new, method = c("scores", "lambda")), "must be a single string")
  expect_error(predict(sfd, new[, "x", drop = FALSE]), "must have the coordinates of the object")

  # En los errores sobre el modelo se muestra la tabla resumen
  add_model(sfd, toSpatFD(fitted["a.PC1"]), name = "incomplete")
  expect_error(predict(sfd, new), "several models.*single.*incomplete")
  expect_error(predict(sfd, new, model = "other"), "must be the name of one of the models.*single")
  expect_error(predict(sfd, new, model = "incomplete", variable = "a"),
               "no variogram for the scores: a.PC2.*single")

  # El kriging simple necesita modelos con meseta
  add_model(sfd, toSpatFD(list(a.PC1 = gstat::vgm(2, "Lin", 0), a.PC2 = fitted$a.PC2)),
                name = "unbounded")
  expect_error(predict(sfd, new, model = "unbounded"), "has no sill")
})

test_that("add_model() registers the models in place", {
  sfd <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")
  expect_equal(nrow(sfd$models), 0)

  # Nombres secuenciales por defecto, o el que dé el usuario
  add_model(sfd, toSpatFD(fitted))
  add_model(sfd, toSpatFD(fitted["a.PC1"]), name = "only PC1")
  add_model(sfd, toSpatFD(fitted))

  expect_equal(rownames(sfd$models), c("model 1", "only PC1", "model 3"))
  expect_named(sfd$models, c("variables", "structures", "no. scores", "cross", "sill"))
  expect_equal(sfd$models$structures, c("nugget, exponential, spherical", "nugget, exponential",
                                        "nugget, exponential, spherical"))
  expect_equal(sfd$models[["no. scores"]], c(2, 1, 2))
  expect_equal(sfd$models$sill, c(33, 32, 33))
  expect_equal(sfd[["models"]], sfd$models)

  # Los modelos compilados quedan en el atributo
  compiled <- attr(sfd, "models")[["model 1"]]
  expect_named(compiled$variograms, c("a.PC1", "a.PC2"))
  expect_equal(compiled$variograms$a.PC1(c(0, 1e6)), c(0, 32))

  # Un nombre repetido reemplaza el modelo sin cambiar su posición
  expect_message(add_model(sfd, toSpatFD(fitted["a.PC2"]), name = "only PC1"), "Replacing")
  expect_equal(rownames(sfd$models), c("model 1", "only PC1", "model 3"))
  expect_equal(sfd$models$structures[2], "spherical")
})

test_that("add_model() translates fitted gstat models on the way in", {
  sfd <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")

  add_model(sfd, fitted, name = "gstat list")
  add_model(sfd, toSpatFD(fitted), name = "table")

  expect_equal(attr(sfd, "models")[["gstat list"]]$table, attr(sfd, "models")[["table"]]$table)
  expect_equal(predict(sfd, new, model = "gstat list")$a$pred$coefs,
               predict(sfd, new, model = "table")$a$pred$coefs)

  expect_error(add_model(sfd, 1), "No SpatFD translator")
})

test_that("add_model() validates the models", {
  sfd <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")

  expect_error(add_model(toSpatFD(fitted)), "only available for SpatFD objects")
  expect_error(add_model(sfd, toSpatFD(list(z.PC1 = fitted$a.PC1))), "do not match the scores")

  df <- toSpatFD(fitted)
  df$model[1] <- "other"
  expect_error(add_model(sfd, df), "Unknown variogram models: other")
  expect_equal(nrow(sfd$models), 0)
})

test_that("models are kept when objects are combined", {
  x <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")
  y <- spatfd(sim_data(sites[1:8], 2), coords, basis = "fourier", nbasis = 11, nharm = 1, name = "b")
  add_model(x, toSpatFD(fitted))
  add_model(y, toSpatFD(list(b.PC1 = gstat::vgm(20, "Exp", 0.5))))

  expect_message(xy <- x + y, "renamed to 'model 2'")
  expect_equal(rownames(xy$models), c("model 1", "model 2"))
  expect_equal(xy$models$variables, c("a", "b"))
  expect_false(any(xy$models$cross))

  # El objeto combinado tiene su propio registro
  add_model(xy, toSpatFD(fitted))
  expect_equal(nrow(xy$models), 3)
  expect_equal(nrow(x$models), 1)
})

test_that("the nugget variogram is zero at the origin", {
  expect_equal(models$nugget(c(0, 1, 5), sigma = 2), c(0, 2, 2))
})
