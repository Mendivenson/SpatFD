h <- seq(10, 1000, length.out = 40)

# Evalúa un modelo de gstat pasando por toSpatFD() + la compilación
compiled <- function(m) compile_model(toSpatFD(list(v = m)))$v
reference <- function(m) gstat::variogramLine(m, dist_vector = h)$gamma

test_that("compiled models reproduce the gstat variogram", {
  gstat_models <- list(
    Nug = gstat::vgm(2, "Nug", 0),
    Exp = gstat::vgm(2, "Exp", 300),
    Sph = gstat::vgm(2, "Sph", 300),
    Gau = gstat::vgm(2, "Gau", 300),
    Exc = gstat::vgm(2, "Exc", 300, kappa = 1.3),
    Mat = gstat::vgm(2, "Mat", 300, kappa = 1.3),
    Ste = gstat::vgm(2, "Ste", 300, kappa = 1.3),
    Lin = gstat::vgm(2, "Lin", 0),
    Per = gstat::vgm(2, "Per", 300),
    Wav = gstat::vgm(2, "Wav", 300),
    Hol = gstat::vgm(2, "Hol", 300),
    Pow = gstat::vgm(2, "Pow", 1.5)
  )

  for (code in names(gstat_models)) {
    m <- gstat_models[[code]]
    expect_equal(compiled(m)(h), reference(m), info = code)
  }
})

test_that("nested structures are added up", {
  m <- gstat::vgm(2, "Exp", 300, nugget = 0.5)
  expect_equal(compiled(m)(h), reference(m))
})

test_that("toSpatFD() returns the tabular schema", {
  fitted <- list(
    "X.PC1" = gstat::vgm(2, "Exp", 300, nugget = 0.5),
    "X.PC2" = gstat::vgm(1, "Sph", 150)
  )
  df <- toSpatFD(fitted)

  expect_named(df, c("var", "model", "rango", "sigma", "alpha1"))
  expect_equal(df$var, c("X.PC1", "X.PC1", "X.PC2"))
  expect_equal(df$model, c("nugget", "exponential", "spherical"))
  expect_named(compile_model(df), c("X.PC1", "X.PC2"))
})

test_that("toSpatFD() rejects unsupported objects", {
  expect_error(toSpatFD(1), "No SpatFD translator")
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
  expect_no_error(models$exponential(h, rango = 3))

  # Solo se revisan los parámetros que el modelo usa
  expect_no_error(models$nugget(h, rango = 0))
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
  df <- data.frame(var = "v", model = "exponential",
                   rango = -300, sigma = 2, alpha1 = NA)
  expect_error(compile_model(df), 'Model "exponential": `rango` must be')
})
