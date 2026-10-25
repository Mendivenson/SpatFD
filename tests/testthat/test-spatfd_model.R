sim_data <- function(sites, seed) {
  set.seed(seed)
  dat <- sapply(sites, function(s) {
    rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(1) * cos(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1)
  })
  colnames(dat) <- sites
  dat
}

sites  <- paste0("s", 1:12)
set.seed(1)
coords <- cbind(x = runif(12), y = runif(12))
rownames(coords) <- sites

a <- suppressMessages(spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 3, name = "a"))
b <- suppressMessages(spatfd(sim_data(sites[1:8], 2), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "b"))

direct <- list(a.PC1 = gstat::vgm(20, "Exp", 0.4, nugget = 1),
               a.PC2 = gstat::vgm(2, "Sph", 0.6))

# ============================ FORMAS DE ENTRADA ===============================
test_that("a list of gstat models builds the model", {
  m <- spatfd_model(a, model = direct)
  expect_s3_class(m, "SpatFDmodel")
  expect_named(m, c("variance", "variograms"))
  expect_equal(m$variograms$variogram, c("a.PC1", "a.PC2"))
  expect_equal(m$variograms$source, c("gstat", "gstat"))
  expect_equal(m$variograms$sill, c(21, 2))
})

test_that("a fitted gstat object builds the model, with cross-variograms", {
  g   <- to_gstat(a, b, pcs = c(a = 1, b = 1))
  fit <- suppressWarnings(gstat::fit.lmc(gstat::variogram(g), g, gstat::vgm(1, "Exp", 0.5, nugget = 0.1), fit.method = 6, correct.diagonal = 1.01))
  m   <- spatfd_model(a, b, model = fit)

  expect_setequal(m$variograms$variogram, c("a.PC1", "b.PC1", "a.PC1.b.PC1"))
  expect_equal(m$variograms$type[m$variograms$variogram == "a.PC1.b.PC1"], "cross")
  expect_equal(unname(m$variance[, "no. PC"]), c(1, 1))
})

test_that("a data frame and a list of functions build the same model as gstat", {
  table <- data.frame(var   = c("a.PC1", "a.PC1", "a.PC2"),
                      model = c("nugget", "exponential", "spherical"),
                      rango = c(NA, 0.4, 0.6),
                      sigma = c(1, 20, 2))
  functions <- list(a.PC1 = function(x) 1 * (x > 0) + 20 * (1 - exp(-x / 0.4)),
                    a.PC2 = function(x) ifelse(x < 0.6, 2 * (1.5 * x / 0.6 - 0.5 * (x / 0.6)^3), 2))

  D  <- as.matrix(dist(coords))
  m1 <- spatfd_model(a, model = direct)
  m2 <- spatfd_model(a, model = table)
  m3 <- spatfd_model(a, model = functions)

  for (id in c("a.PC1", "a.PC2")) {
    expect_equal(pair_variogram(m2, id, id)(D), pair_variogram(m1, id, id)(D))
    expect_equal(pair_variogram(m3, id, id)(D), pair_variogram(m1, id, id)(D))
  }
  expect_equal(m2$variograms$sill, c(21, 2))
  expect_equal(m3$variograms$sill, c(21, 2))        # La meseta es f(Inf)
})

test_that("the sill of a function is searched in cascade", {
  sill <- function(f) spatfd_model(a, model = list(a.PC1 = f))$variograms$sill

  # 1. f(Inf): finito es la meseta, infinito es un modelo sin meseta
  expect_equal(sill(function(x) 2 * (1 - exp(-x / 0.3))), 2)
  expect_true(is.na(sill(function(x) 2 * x)))

  # 2. Límite numérico: sin(Inf) no está definido, pero el modelo wave converge
  wave <- function(x) ifelse(x == 0, 0, 2 * (1 - sin(pi * x / 0.3) / (pi * x / 0.3)))
  expect_message(expect_equal(sill(wave), 2), "limit at large distances")

  # 3. Oscilación acotada: el promedio, aunque no sea simétrica. Un periódico
  # isotrópico no es un modelo válido en el plano, así que se prueba la
  # búsqueda de la meseta directamente
  periodic <- function(x) 2 * (1 - cos(2 * pi * x / 0.3)) + (1 - cos(4 * pi * x / 0.3))
  expect_message(expect_equal(function_sill(periodic, "p", 1), 3, tolerance = 0.01), "oscillates")

  # Si crece sin límite no tiene meseta
  expect_message(expect_true(is.na(function_sill(function(x) x / 0.3 + sin(x / 0.3), "g", 1))),
                 "grows without a limit")
})

test_that("an invalid model is stopped when it is built", {
  periodic <- function(x) 2 * (1 - cos(2 * pi * x / 0.3))
  expect_error(suppressMessages(spatfd_model(a, model = list(a.PC1 = periodic))), "The model is not valid")

  # Un variograma cruzado más fuerte de lo que permiten los directos
  lmc <- list(a.PC1 = gstat::vgm(1, "Exp", 0.5), b.PC1 = gstat::vgm(1, "Exp", 0.5),
              a.PC1.b.PC1 = gstat::vgm(5, "Exp", 0.5))
  expect_error(spatfd_model(a, b, model = lmc), "The model is not valid")

  # La revisión por estructura de un modelo de corregionalización: una meseta
  # parcial cruzada imposible solo en la estructura exponencial
  lmc <- list(a.PC1 = gstat::vgm(1, "Exp", 0.5, nugget = 1), b.PC1 = gstat::vgm(1, "Exp", 0.5, nugget = 1),
              a.PC1.b.PC1 = gstat::vgm(2, "Exp", 0.5, nugget = 0))
  expect_error(spatfd_model(a, b, model = lmc), "The model is not valid")
})

test_that("unsupported inputs are rejected", {
  expect_error(spatfd_model(a), "must be given in 'model'")
  expect_error(spatfd_model(a, model = 1), "'model' must be a fitted gstat object")
  expect_error(spatfd_model(a, model = gstat::vgm(1, "Exp", 1)), "A single variogram model has no id")
  expect_error(spatfd_model(a, model = list(gstat::vgm(1, "Exp", 1))), "'model' must be a fitted gstat object")
  expect_error(spatfd_model(a, model = list(a.PC1 = function(x) x + 1)), "0 at distance 0")
  expect_error(spatfd_model(a, model = list(a.PC1 = function(x, r = 1) x / r)), "a single argument")
})

test_that("anisotropic gstat models are not supported yet", {
  expect_error(spatfd_model(a, model = list(a.PC1 = gstat::vgm(1, "Exp", 1, anis = c(30, 0.5)))),
               "Anisotropic variogram models are not supported yet by SpatFD \\(model of 'a.PC1'\\)")
})

# ================================== IDS =======================================
test_that("the ids must match the scores of the objects", {
  expect_error(spatfd_model(a, model = list(c.PC1 = gstat::vgm(1, "Exp", 1))),
               "do not match the scores of the objects: c.PC1")
  expect_error(spatfd_model(a, b, model = direct), "no variogram for the variable 'b'")
  expect_error(spatfd_model(a, model = list(a.PC2 = gstat::vgm(1, "Exp", 1))), "without gaps. Found: PC2")
  expect_error(spatfd_model(a, model = list(a.PC1 = gstat::vgm(1, "Exp", 1), a.PC1.a.PC2 = gstat::vgm(1, "Exp", 1))),
               "no direct variogram: a.PC2")
  expect_error(spatfd_model(a, a, model = direct), "Several objects are named 'a'")
})

test_that("cross ids are read in any order and with dots in the names", {
  pm <- suppressMessages(spatfd(sim_data(sites, 3), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "PM2.5"))
  m  <- spatfd_model(a, pm, model = list(a.PC1 = gstat::vgm(1, "Exp", 1), PM2.5.PC1 = gstat::vgm(1, "Exp", 1),
                                         PM2.5.PC1.a.PC1 = gstat::vgm(0.5, "Exp", 1)))
  expect_equal(m$variograms$type, c("direct", "direct", "cross"))
  expect_false(is.null(pair_variogram(m, "a.PC1", "PM2.5.PC1")))
  expect_null(pair_variogram(m, "a.PC1", "a.PC2"))
})

# ============================ COPIAS TRUNCADAS ================================
test_that("the model keeps each variable with its first K components", {
  m    <- spatfd_model(a, model = direct[1])
  copy <- get_model_data(m)$a

  expect_equal(colnames(copy$scores), c("x", "y", "a.PC1"))
  expect_equal(ncol(get_pca(copy)$scores), 1)
  expect_equal(ncol(get_pca(copy)$harmonics$coefs), 1)
  expect_equal(m$variance$variance, get_pca(a)$varprop[1])
  expect_equal(attr(copy, "name"), "a")
  expect_equal(attr(copy, "internal")$data, attr(a, "internal")$data)
})

# ================================ MÉTODOS =====================================
test_that("print() and summary() show the tables", {
  m <- spatfd_model(a, model = direct)
  expect_output(print(m), "Variograms")
  expect_output(print(summary(m)), "Variables")
})

test_that("plot() draws the matrix of variograms", {
  g   <- to_gstat(a, b, pcs = c(a = 1, b = 1))
  fit <- suppressWarnings(gstat::fit.lmc(gstat::variogram(g), g, gstat::vgm(1, "Exp", 0.5, nugget = 0.1), fit.method = 6, correct.diagonal = 1.01))
  m   <- spatfd_model(a, b, model = fit)

  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(m))
  expect_invisible(plot(spatfd_model(a, model = direct), cutoff = 0.5))
})

test_that("plot() draws several models, with all their scores", {
  g   <- to_gstat(a, b, pcs = c(a = 1, b = 1))
  fit <- suppressWarnings(gstat::fit.lmc(gstat::variogram(g), g, gstat::vgm(1, "Exp", 0.5, nugget = 0.1), fit.method = 6, correct.diagonal = 1.01))
  m1  <- spatfd_model(a, b, model = fit)
  m2  <- spatfd_model(a, model = direct)              # Otra K para 'a' y sin 'b'

  # Todas las parejas de scores de los dos modelos
  data <- models_data(list(m1, m2))
  expect_equal(names(data), c("a", "b"))
  expect_equal(ncol(get_pca(data$a)$scores), 2)

  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(m1, m2, cutoff = 0.8))
  expect_invisible(plot(m1, second = m2))
  expect_error(plot(m1, 0.5), "unnamed arguments must be SpatFDmodel")

  # Variables con el mismo nombre y scores distintos
  other <- suppressMessages(spatfd(sim_data(sites, 9), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a"))
  m3    <- spatfd_model(other, model = direct)
  expect_warning(plot(m2, m3), "different scores in the models")
})
