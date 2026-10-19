sim_data <- function(sites, seed) {
  set.seed(seed)
  t <- 1:50
  dat <- sapply(sites, function(i) rnorm(1, 5) * sin(2 * pi * t / 50) + rnorm(50, sd = 0.1))
  colnames(dat) <- sites
  dat
}

sites  <- paste0("s", 1:10)
set.seed(1)
coords <- cbind(x = runif(10), y = runif(10))
rownames(coords) <- sites

a <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 3, name = "a")
a <- add_model(a, list(a.PC1 = gstat::vgm(30, "Exp", 0.5),
                       a.PC2 = gstat::vgm(1, "Exp", 0.5),
                       a.PC3 = gstat::vgm(0.1, "Exp", 0.5)))
pred <- predict(a, cbind(x = c(0.2, 0.7), y = c(0.5, 0.5)))

# Cuenta las páginas nuevas que se abren con el dispositivo en modo 'ask',
# es decir, las veces que R pediría <Enter> en una sesión interactiva
pages_asking <- function(code) {
  asked <- 0
  setHook("before.plot.new", function() if (par("ask")) asked <<- asked + 1)
  on.exit(setHook("before.plot.new", NULL, "replace"))
  force(code)
  asked
}

test_that("no plot asks for <Enter>, even with the device in 'ask' mode", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  par(ask = TRUE)

  expect_equal(pages_asking(suppressMessages(plot(a))), 0)
  expect_equal(pages_asking(suppressMessages(plot(a, mode = "pca"))), 0)
  expect_equal(pages_asking(plot(pred)), 0)
  expect_equal(pages_asking(availablemodels()), 0)

  if (requireNamespace("sf", quietly = TRUE)) {
    square <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))))
    expect_equal(pages_asking(plot(predict(a, square, n = 100))), 0)
  }

  # El modo del dispositivo se respeta: se restaura al terminar
  expect_true(par("ask"))
})

test_that("the principal components are drawn on a single page", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  pages <- 0
  setHook("before.plot.new", function() pages <<- pages + 1)
  on.exit(setHook("before.plot.new", NULL, "replace"), add = TRUE)

  plot(a, mode = "pca")
  expect_equal(pages, 3)                      # tres paneles, una sola página
  expect_equal(par("mfrow"), c(1, 1))         # la división de la página se restaura
})

test_that("plot() returns the object invisibly, without printing", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  expect_invisible(plot(a))
  expect_silent(plot(a, variable = "a"))
})
