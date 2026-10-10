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

a <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")
b <- spatfd(sim_data(sites[1:8], 2), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "b")

test_that("spatfd() builds the expected object", {
  expect_s3_class(a, "SpatFD")
  expect_equal(var_names(a), "a")
  expect_equal(colnames(a$scores), c("x", "y", "a.PC1", "a.PC2"))
  expect_equal(rownames(a$scores), sites)
  expect_equal(rownames(a$variance), "a")
})

test_that("spatfd() validates its arguments", {
  expect_error(spatfd(as.data.frame(sim_data(sites, 1)), coords), "'dat' must be a non-empty matrix")
  expect_error(spatfd(sim_data(sites, 1), coords, basis = "other", nbasis = 11, nharm = 2),
               "'basis' must be either")
})

test_that("`+` joins variables and sites", {
  ab <- a + b

  expect_s3_class(ab, "SpatFD")
  expect_equal(var_names(ab), c("a", "b"))
  expect_equal(nrow(ab$scores), 10)
  expect_true(all(is.na(ab$scores[c("s9", "s10"), "b.PC1"])))
  expect_equal(rownames(ab$variance), c("a", "b"))
})

test_that("togstat() exports the selected scores", {
  ab <- a + b

  g <- togstat(ab)
  expect_s3_class(g, "gstat")
  expect_named(g$data, c("a.PC1", "a.PC2", "b.PC1", "b.PC2"))
  expect_equal(nrow(g$data$b.PC1$data), 8)

  expect_named(togstat(ab, vars = "a")$data, c("a.PC1", "a.PC2"))
  expect_error(togstat(ab, vars = "z"), "Variables not found")
})

test_that("combining objects does not mix up sites with similar names", {
  # 's1' solo está en la segunda variable: no debe tomar los scores de 's10'
  x  <- spatfd(sim_data(sites[-1], 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "x")
  y  <- spatfd(sim_data(sites, 2), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "y")
  xy <- x + y

  expect_true(all(is.na(xy$scores["s1", c("x.PC1", "x.PC2")])))
  expect_equal(xy$scores["s10", c("x.PC1", "x.PC2")], x$scores["s10", c("x.PC1", "x.PC2")])
  expect_equal(xy$scores[sites, c("y.PC1", "y.PC2")], y$scores[sites, c("y.PC1", "y.PC2")])
})
