# Dos variables con estaciones distintas, como PM10 y NO2 en México
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

test_that("to_gstat() exports one entry per score", {
  g <- to_gstat(a, b)
  expect_s3_class(g, "gstat")
  expect_equal(names(g$data), c("a.PC1", "a.PC2", "a.PC3", "b.PC1", "b.PC2"))
})

test_that("each variable keeps its own sites", {
  g <- to_gstat(a, b)
  expect_equal(nrow(g$data$a.PC1$data), 12)
  expect_equal(nrow(g$data$b.PC1$data), 8)
  expect_false(anyNA(g$data$b.PC1$data$b.PC1))
  expect_equal(unname(g$data$a.PC2$data$a.PC2), unname(a$scores[, "a.PC2"]))
})

test_that("pcs selects the first components of each variable", {
  g <- to_gstat(a, b, pcs = c(a = 2))
  expect_equal(names(g$data), c("a.PC1", "a.PC2", "b.PC1", "b.PC2"))

  g <- to_gstat(b, a, pcs = c(a = 1, b = 1))
  expect_equal(names(g$data), c("b.PC1", "a.PC1"))
})

test_that("to_gstat() validates its arguments", {
  expect_error(to_gstat(a, list()), "must be of class SpatFD")
  expect_error(to_gstat(a, a), "Several objects are named 'a'")
  expect_error(to_gstat(a, pcs = c(c = 1)), "not given in '...': c")
  expect_error(to_gstat(a, pcs = c(a = 4)), "Invalid for: a \\(3 components\\)")
  expect_error(to_gstat(a, pcs = 1), "named numeric vector")

  other <- coords
  colnames(other) <- c("lon", "lat")
  c2 <- suppressMessages(spatfd(sim_data(sites, 3), other, basis = "fourier", nbasis = 11, nharm = 2, name = "c"))
  expect_error(to_gstat(a, c2), "same column names")
})

test_that("the export can be fitted by gstat", {
  g <- to_gstat(a, b, pcs = c(a = 1, b = 1))
  v <- gstat::variogram(g)
  expect_setequal(unique(as.character(v$id)), c("a.PC1", "b.PC1", "a.PC1.b.PC1"))
})
