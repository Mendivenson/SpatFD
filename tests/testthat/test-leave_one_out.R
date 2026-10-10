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

new_a <- function(dat = sim_data(sites, 1)) {
  spatfd(dat, coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")
}
new_b <- function() {
  spatfd(sim_data(sites[1:8], 2), coords, basis = "fourier", nbasis = 11, nharm = 1, name = "b")
}

exponential <- list(a.PC1 = gstat::vgm(30, "Exp", 0.5, nugget = 2), a.PC2 = gstat::vgm(1, "Exp", 0.5))
spherical   <- list(a.PC1 = gstat::vgm(30, "Sph", 0.8, nugget = 2), a.PC2 = gstat::vgm(1, "Sph", 0.8))

test_that("leave_one_out() compares every model with every method by default", {
  a <- new_a()
  add_model(a, exponential, name = "exponential")
  add_model(a, spherical, name = "spherical")

  loo <- leave_one_out(a)

  expect_s3_class(loo, "SpatFDloo")
  expect_named(loo, c("summary", "errors", "pred"))
  expect_named(loo$summary, c("model", "method", "crossed", "mean", "sd", "max"))
  expect_equal(nrow(loo$summary), 4)
  expect_setequal(paste(loo$summary$model, loo$summary$method),
                  c("exponential scores", "exponential lambda", "spherical scores", "spherical lambda"))
  expect_false(is.unsorted(loo$summary$mean))                 # del menor al mayor error
  expect_false(any(loo$summary$crossed))

  expect_equal(dim(loo$errors), c(10, 4))
  expect_equal(rownames(loo$errors), sites)
  expect_true(all(loo$errors > 0))
  expect_equal(loo$summary$mean, unname(sort(colMeans(loo$errors))))

  expect_s3_class(loo$pred[["exponential / scores"]], "fd")
  expect_equal(colnames(loo$pred[["exponential / scores"]]$coefs), sites)

  # El registro de modelos del objeto original no cambia
  expect_equal(rownames(a$models), c("exponential", "spherical"))
})

test_that("each error is the one of predicting the site after repeating the analysis without it", {
  dat <- sim_data(sites, 1)
  a   <- new_a(dat)
  add_model(a, exponential, name = "exponential")
  loo <- leave_one_out(a)

  without <- new_a(dat[, sites != "s3"])
  add_model(without, exponential, name = "exponential")

  basis <- attr(a, "fd")$a$basis
  gram  <- fda::eval.penalty(basis, fda::int2Lfd(0))

  for (method in c("scores", "lambda")) {
    pred <- predict(without, coords["s3", , drop = FALSE], method = method)$a$pred
    difference <- attr(a, "fd")$a$coefs[, "s3"] - pred$coefs[, 1]

    label <- paste("exponential /", method)
    expect_equal(loo$errors["s3", label], sqrt(drop(crossprod(difference, gram %*% difference))),
                 info = method)
    expect_equal(unname(loo$pred[[label]]$coefs[, "s3"]), unname(pred$coefs[, 1]), info = method)
  }
})

test_that("the design chooses the combinations to compare", {
  a <- new_a()
  add_model(a, exponential, name = "exponential")
  add_model(a, spherical, name = "spherical")

  design <- data.frame(model = c("spherical", "exponential"), method = c("scores", "lambda"))
  loo    <- leave_one_out(a, design = design)

  expect_equal(nrow(loo$summary), 2)
  expect_setequal(colnames(loo$errors), c("spherical / scores", "exponential / lambda"))

  expect_error(leave_one_out(a, design = data.frame(model = "other", method = "scores")),
               "Models not found in the object: other.*exponential")
  expect_error(leave_one_out(a, design = data.frame(model = "spherical", method = "other")),
               "'method' must be")
  expect_error(leave_one_out(a, design = data.frame(model = "spherical")), "'design' must be")
  expect_error(leave_one_out(a, variable = "z"), "'variable' must be")
  expect_error(leave_one_out(new_a()), "no combinations to compare")
})

test_that("multivariate models are cross-validated one variable at a time", {
  ab <- new_a() + new_b()

  direct  <- c(exponential, list(b.PC1 = gstat::vgm(20, "Exp", 0.5, nugget = 2)))
  crossed <- toSpatFD(direct)
  crossed <- rbind(crossed,
                   transform(crossed[crossed$var == "b.PC1", ], var = "a.PC1.b.PC1", sigma = c(0.5, -10)))
  add_model(ab, direct, name = "direct")
  add_model(ab, crossed, name = "crossed")

  expect_message(loo <- leave_one_out(ab), "Cross-validating the variable a")

  # lambda no admite modelos con variogramas cruzados: no entra por defecto
  expect_setequal(paste(loo$summary$model, loo$summary$method),
                  c("direct scores", "direct lambda", "crossed scores"))
  expect_equal(loo$summary$crossed[loo$summary$model == "crossed"], TRUE)
  expect_equal(rownames(loo$errors), sites)

  # La otra variable se valida sobre sus propios sitios
  loo_b <- leave_one_out(ab, variable = "b")
  expect_equal(rownames(loo_b$errors), sites[1:8])

  expect_error(leave_one_out(ab, design = data.frame(model = "crossed", method = "lambda")),
               "has cross-variograms")
})
