skip_if_not_installed("plotly")
skip_if_not_installed("sf")

sim_data <- function(sites, time) {
  set.seed(1)
  dat <- data.frame(time = time)
  for (s in sites) dat[[s]] <- rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1)
  dat
}

sites  <- paste0("s", 1:10)
set.seed(1)
coords <- cbind(x = runif(10), y = runif(10))
rownames(coords) <- sites
square <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))))
table  <- data.frame(var = "a.PC1", model = "exponential", rango = 0.5, sigma = 30)

a <- spatfd(sim_data(sites, 1:50), coords, basis = "fourier", nbasis = 11, nharm = 1, name = "a")
p <- predict(spatfd_model(a, model = table), square, grid_size = 100)

test_that("plot_animated() makes one frame per time", {
  fig <- plot_animated(p, times = 5)
  expect_s3_class(fig, "plotly")
  expect_length(fig$x$frames, 5)
  expect_equal(vapply(fig$x$frames, `[[`, character(1), "name"), c("1", "13.25", "25.5", "37.75", "50"))

  # Solo se anima el mapa de la predicción
  expect_equal(fig$x$frames[[1]]$traces, list(0))
})

test_that("the times of the slider can be given in the format of the data", {
  expect_equal(vapply(plot_animated(p, times = c(10, 3))$x$frames, `[[`, character(1), "name"), c("3", "10"))

  daily <- spatfd(sim_data(sites, as.Date("2026-01-01") + 0:49), coords, basis = "fourier", nbasis = 11,
                  nharm = 1, name = "a")
  pd    <- predict(spatfd_model(daily, model = table), square, grid_size = 100)
  fig   <- plot_animated(pd, times = as.Date(c("2026-01-05", "2026-02-01")))
  expect_equal(vapply(fig$x$frames, `[[`, character(1), "name"), c("2026-01-05", "2026-02-01"))
})

test_that("plot_animated() validates its arguments", {
  expect_error(plot_animated(p, times = 1), "whole number of at least 2")
  expect_error(plot_animated(p, times = c(0, 10)), "within the times of the curves")
  expect_error(plot_animated(p, variable = "b"), "Variable b was not predicted")

  sites_pred <- predict(spatfd_model(a, model = table), cbind(x = 0.5, y = 0.5))
  expect_error(plot_animated(sites_pred), "only available for predictions made over a map")
})
