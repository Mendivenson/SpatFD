h <- c(0, seq(10, 1000, length.out = 40))

# Nombres de gstat de los modelos del catálogo que lo imitan
gstat_names <- c(Nug = "nugget", Exp = "exponential", Sph = "spherical", Gau = "gaussian", Exc = "stable",
                 Mat = "matern", Ste = "stein's matern", Cir = "circular", Lin = "linear", Bes = "bessel",
                 Pen = "pentamodel", Per = "periodic", Wav = "wave", Hol = "hole", Log = "logarithmic",
                 Pow = "power", Spl = "spline", Leg = "legendre")

# Un modelo de gstat por las dos rutas: directo con gstat y como data frame del catálogo
by_gstat <- function(m) from_gstat(list(v = m))$v
by_table <- function(m) {
  from_table(data.frame(var = "v", model = unname(gstat_names[as.character(m$model)]),
                        rango = m$range, sigma = m$psill, alpha1 = m$kappa))$v
}

gstat_models <- list(
  Nug = gstat::vgm(2, "Nug", 0),
  Exp = gstat::vgm(2, "Exp", 300),
  Sph = gstat::vgm(2, "Sph", 300),
  Gau = gstat::vgm(2, "Gau", 300),
  Exc = gstat::vgm(2, "Exc", 300, kappa = 1.3),
  Mat = gstat::vgm(2, "Mat", 300, kappa = 1.3),
  Ste = gstat::vgm(2, "Ste", 300, kappa = 1.3),
  Cir = gstat::vgm(2, "Cir", 300),
  Lin = gstat::vgm(2, "Lin", 0),
  Lin_with_range = gstat::vgm(2, "Lin", 300),
  Bes = gstat::vgm(2, "Bes", 300),
  Pen = gstat::vgm(2, "Pen", 300),
  Per = gstat::vgm(2, "Per", 300),
  Wav = gstat::vgm(2, "Wav", 300),
  Hol = gstat::vgm(2, "Hol", 300),
  Log = gstat::vgm(2, "Log", 300),
  Pow = gstat::vgm(2, "Pow", 1.5),
  Spl = gstat::vgm(2, "Spl", 300),
  Leg = gstat::vgm(2, "Leg", 0.5),
  Nested = gstat::vgm(2, "Exp", 300, nugget = 0.5)
)

test_that("the catalog reproduces the gstat variograms", {
  for (code in names(gstat_models)) {
    m <- gstat_models[[code]]
    expect_equal(by_table(m)$variogram(h), by_gstat(m)$variogram(h), info = code)
  }
})

test_that("both routes agree on which models have a sill", {
  for (code in names(gstat_models)) {
    m <- gstat_models[[code]]
    expect_equal(by_table(m)$sill, by_gstat(m)$sill, info = code)
  }
  expect_true(is.na(by_gstat(gstat_models$Lin)$sill))
  expect_true(is.na(by_gstat(gstat_models$Pow)$sill))
  expect_equal(by_gstat(gstat_models$Lin_with_range)$sill, 2)
  expect_equal(by_gstat(gstat_models$Nested)$sill, 2.5)
})

test_that("the compiled variograms keep the shape of the distances", {
  D <- matrix(c(0, 100, 100, 0), 2)
  for (route in list(by_gstat, by_table)) {
    out <- route(gstat_models$Nested)$variogram(D)
    expect_equal(dim(out), c(2, 2))
    expect_equal(out[1, 1], 0)
  }
  expect_length(by_gstat(gstat_models$Exp)$variogram(numeric(0)), 0)
})

test_that("models reject parameters outside their domain", {
  expect_error(models$exponential(h, rango = 0), "`rango` must be")
  expect_error(models$exponential(h, rango = -300), "`rango` must be")
  expect_error(models$exponential(h, rango = c(1, 2)), "`rango` must be")
  expect_error(models$exponential(h, sigma = NA_real_), "`sigma` must be")
  expect_error(models$exponential(h, sigma = "2"), "`sigma` must be")
  expect_error(models$matern(h, alpha1 = 0), "`alpha1` must be")
  expect_error(models$exponential(-h), "`h` must be")
  expect_error(models$exponential(c(h, NA)), "`h` must be")

  # Cotas propias de cada modelo
  expect_error(models$stable(h, alpha1 = 2.5), "in (0, 2]", fixed = TRUE)
  expect_error(models$power(h, rango = 3), "in (0, 2]", fixed = TRUE)
  expect_no_error(models$stable(h, alpha1 = 2))

  # El modelo lineal admite rango cero (sin meseta), pero no negativo
  expect_no_error(models$linear(h, rango = 0))
  expect_error(models$linear(h, rango = -1), "`rango` must be")

  # Mesetas parciales negativas: válidas en variogramas cruzados
  expect_no_error(models$exponential(h, sigma = -2))
})

test_that("every model accepts its default parameters", {
  for (model in names(models)) {
    expect_no_error(models[[model]](h))
    expect_length(models[[model]](numeric(0)), 0)
  }
})

test_that("invalid parameters are reported when compiling", {
  df <- data.frame(var = "v", model = "exponential", rango = -300, sigma = 2)
  expect_error(from_table(df), 'Model "exponential": `rango` must be')
  expect_error(from_table(data.frame(var = "v", model = "other", sigma = 2)), "Unknown variogram models: other")
  expect_error(from_table(data.frame(var = "v", sigma = 2)), "needs the columns")
})

test_that("availablemodels() lists and plots the models", {
  expect_output(availablemodels(txt = TRUE), "exponential")
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(availablemodels())
})
