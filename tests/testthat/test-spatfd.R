# Datos simulados: una columna 'time' y una columna por sitio
sim_data <- function(sites, seed, time = 1:50) {
  set.seed(seed)
  dat <- data.frame(time = time)
  for (s in sites) {
    dat[[s]] <- rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(1) * cos(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1)
  }
  dat
}

sites  <- paste0("s", 1:10)
set.seed(1)
coords <- cbind(x = runif(10), y = runif(10))
rownames(coords) <- sites

dat <- sim_data(sites, 1)
a   <- spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")

# ================================ OBJETO ======================================
test_that("spatfd() builds the expected object", {
  expect_s3_class(a, "SpatFD")
  expect_named(a, c("scores", "variance", "missing"))
  expect_equal(attr(a, "name"), "a")
  expect_equal(colnames(a$scores), c("x", "y", "a.PC1", "a.PC2"))
  expect_equal(rownames(a$scores), sites)
  expect_equal(a$scores[, c("x", "y")], coords)
  expect_equal(rownames(a$variance), "a")
  expect_equal(unname(a$missing), rep(0, 10))
})

test_that("the default name is the expression given as 'dat'", {
  expect_equal(attr(spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2), "name"), "dat")
})

test_that("spatfd() validates its arguments", {
  expect_error(spatfd(list(1), coords), "'dat' must be a non-empty data frame or matrix")
  expect_error(spatfd(dat, coords, basis = "other", nbasis = 11, nharm = 2), "'basis' must be either")
  expect_error(spatfd(dat, coords[1:5, ], basis = "fourier", nbasis = 11, nharm = 2), "without coordinates: s6")
  expect_error(spatfd(transform(dat, s1 = "a"), coords, basis = "fourier", nbasis = 11, nharm = 2),
               "Not numeric: s1")
  expect_error(spatfd(transform(dat, time = 1), coords, basis = "fourier", nbasis = 11, nharm = 2),
               "The times must be unique")

  same <- coords
  same["s7", ] <- same["s3", ]
  expect_error(spatfd(dat, same, basis = "fourier", nbasis = 11, nharm = 2),
               "Sites with the same coordinates: s3 and s7")
})

# =========================== DATA FRAME Y MATRIZ ==============================
test_that("a matrix gives the same object as a data frame", {
  expect_equal(spatfd(as.matrix(dat), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a"), a)
})

test_that("coordinates can be a data frame with a column of sites", {
  cdf <- data.frame(station = sites, x = coords[, "x"], y = coords[, "y"])
  b   <- spatfd(dat, cdf, id = "station", basis = "fourier", nbasis = 11, nharm = 2, name = "a")
  expect_equal(b$scores, a$scores)
  expect_error(spatfd(dat, cdf, id = "other", basis = "fourier", nbasis = 11, nharm = 2), "No column 'other'")
})

# ================================ TIEMPOS =====================================
test_that("without a column of times the rows are equally spaced", {
  expect_message(b <- spatfd(dat[, -1], coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a"),
                 "equally spaced times 1, ..., 50")
  expect_equal(b$scores, a$scores)
})

test_that("the column of times can have another name", {
  other <- dat
  names(other)[1] <- "hour"
  expect_equal(spatfd(other, coords, time = "hour", basis = "fourier", nbasis = 11, nharm = 2, name = "a"), a)
})

test_that("the rows are sorted by time", {
  o <- sample(50)
  b <- spatfd(dat[o, ], coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")
  expect_equal(b$scores, a$scores)
  expect_equal(attr(b, "internal")$data, attr(a, "internal")$data)
})

test_that("dates become days, with the time of day as a fraction", {
  expect_equal(as_days(as.Date("1970-01-11")), 10)
  expect_equal(as_days(as.POSIXct("1970-01-11 06:00:30", tz = "UTC")), 10.25 + 30 / 86400)

  # Siempre en UTC: una fecha es las 00:00 UTC de ese día, y Bogotá va 5 horas atrás
  expect_equal(as_days(as.POSIXct("2026-10-25 00:00", tz = "UTC")), as_days(as.Date("2026-10-25")))
  expect_equal(as_days(as.POSIXct("2026-10-25 00:00", tz = "America/Bogota")), as_days(as.Date("2026-10-25")) + 5 / 24)

  # El cambio de horario no repite horas: 01:30 ocurre dos veces en Nueva York el 2026-11-01
  ny <- as.POSIXct("2026-11-01 00:30", tz = "America/New_York") + 3600 * 0:3
  expect_equal(diff(as_days(ny)), rep(1 / 24, 3))

  daily <- sim_data(sites, 1, time = as.Date("2026-01-01") + 0:49)
  b     <- spatfd(daily, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")
  expect_equal(attr(b, "internal")$time, daily$time)
  expect_equal(attr(b, "internal")$t, as.numeric(daily$time))
})

# =========================== DATOS FALTANTES ==================================
test_that("a site with missing values is smoothed with its observations", {
  miss <- dat
  miss$s2[c(3, 10, 11, 40)] <- NA
  b <- spatfd(miss, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")

  coefs <- attr(b, "internal")$fd$coefs
  expect_false(anyNA(coefs))
  expect_equal(unname(b$missing["s2"]), 8)

  # La curva de s2 es la misma que si se suavizara sola sin esos tiempos
  fdpar <- fda::fdPar(fda::create.fourier.basis(c(1, 50), 11), fda::int2Lfd(2), 0.0002)
  alone <- fda::smooth.basis(dat$time[-c(3, 10, 11, 40)], dat$s2[-c(3, 10, 11, 40)], fdpar)$fd
  expect_equal(unname(coefs[, "s2"]), as.vector(alone$coefs))

  # Las demás curvas no cambian
  expect_equal(coefs[, colnames(coefs) != "s2"], attr(a, "internal")$fd$coefs[, colnames(coefs) != "s2"])
})

test_that("sites with more than 10% of missing values raise a warning", {
  miss <- dat
  miss$s3[1:6] <- NA
  expect_warning(b <- spatfd(miss, coords, basis = "fourier", nbasis = 11, nharm = 2),
                 "more than 10% of missing values: s3 \\(12%\\)")
  expect_true("s3" %in% rownames(b$scores))
})

test_that("only the sites without observations are removed", {
  miss <- dat
  miss$s4 <- NA_real_
  expect_message(b <- spatfd(miss, coords, basis = "fourier", nbasis = 11, nharm = 2),
                 "Sites without observations, removed: s4")
  expect_equal(rownames(b$scores), setdiff(sites, "s4"))
  expect_false("s4" %in% names(b$missing))
})

test_that("rows without observations are dropped", {
  miss <- dat
  miss[c(1, 50), -1] <- NA
  b <- spatfd(miss, coords, basis = "fourier", nbasis = 11, nharm = 2)
  expect_equal(attr(b, "internal")$t, 2:49)
  expect_equal(unname(b$missing), rep(0, 10))
})

test_that("without penalty every site needs at least nbasis observations", {
  miss <- dat
  miss$s5[1:45] <- NA
  expect_error(suppressWarnings(spatfd(miss, coords, basis = "fourier", nbasis = 11, lambda = 0, nharm = 2)),
               "every site needs at least 'nbasis' observations")
})

# ================================ MÉTODOS =====================================
test_that("summary() returns the variance and the missing values", {
  s <- summary(a)
  expect_s3_class(s, "summary.SpatFD")
  expect_equal(s$variance, a$variance)
  expect_output(print(s), "PCA information")
  expect_output(print(a), "a.PC1")
})

test_that("plot() draws both modes", {
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(a))
  expect_invisible(plot(a, mode = "pca"))
  expect_error(plot(a, mode = "other"), "'mode' must be")
})
