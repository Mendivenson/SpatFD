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
  a <- new_a() |>
    add_model(exponential, name = "exponential") |>
    add_model(spherical, name = "spherical")

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

  # Las curvas predichas son matrices tiempos x sitios, como los datos de entrada
  expect_true(is.matrix(loo$pred[["exponential / scores"]]))
  expect_equal(dim(loo$pred[["exponential / scores"]]), c(50, 10))
  expect_equal(colnames(loo$pred[["exponential / scores"]]), sites)

  # El registro de modelos del objeto original no cambia
  expect_equal(rownames(a$models), c("exponential", "spherical"))
})

test_that("each error is the one of predicting the site after repeating the analysis without it", {
  dat <- sim_data(sites, 1)
  a   <- add_model(new_a(dat), exponential, name = "exponential")
  loo <- leave_one_out(a)

  without <- add_model(new_a(dat[, sites != "s3"]), exponential, name = "exponential")

  observed <- get_fd(a, "a")
  gram     <- fda::eval.penalty(observed$basis, fda::int2Lfd(0))

  for (method in c("scores", "lambda")) {
    pred       <- predict(without, coords["s3", , drop = FALSE], method = method)
    difference <- observed$coefs[, "s3"] - predicted_fd(pred, "a")$coefs[, 1]

    label <- paste("exponential /", method)
    expect_equal(loo$errors["s3", label], sqrt(drop(crossprod(difference, gram %*% difference))),
                 info = method)
    expect_equal(unname(loo$pred[[label]][, "s3"]), unname(pred$a$pred[, 1]), info = method)
  }
})

test_that("the design chooses the combinations to compare", {
  a <- new_a() |>
    add_model(exponential, name = "exponential") |>
    add_model(spherical, name = "spherical")

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
  ab <- ab |>
    add_model(direct, name = "direct") |>
    add_model(crossed, name = "crossed")

  expect_message(loo <- leave_one_out(ab), "Cross-validating the variable a")

  # Cada modelo con cada método, también lambda con variogramas cruzados
  expect_setequal(paste(loo$summary$model, loo$summary$method),
                  c("direct scores", "direct lambda", "crossed scores", "crossed lambda"))
  expect_equal(loo$summary$crossed, loo$summary$model == "crossed")
  expect_equal(rownames(loo$errors), sites)

  # La otra variable se valida sobre sus propios sitios
  loo_b <- leave_one_out(ab, variable = "b")
  expect_equal(rownames(loo_b$errors), sites[1:8])
})

test_that("with different time intervals lambda is cross-validated over the common one", {
  short <- spatfd(sim_data(sites[1:8], 2)[1:40, ], coords, t = 1:40,
                  basis = "fourier", nbasis = 11, nharm = 1, name = "b")
  crossed <- toSpatFD(c(exponential, list(b.PC1 = gstat::vgm(20, "Exp", 0.5, nugget = 2))))
  crossed <- rbind(crossed,
                   transform(crossed[crossed$var == "b.PC1", ], var = "a.PC1.b.PC1", sigma = c(0.5, -10)))
  ab <- add_model(new_a() + short, crossed, name = "crossed")

  loo <- leave_one_out(ab, variable = "a")

  # Las curvas de lambda solo cubren el intervalo común; las de scores, todo
  expect_true(all(is.na(loo$pred[["crossed / lambda"]][as.character(41:50), ])))
  expect_false(anyNA(loo$pred[["crossed / lambda"]][as.character(1:40), ]))
  expect_false(anyNA(loo$pred[["crossed / scores"]]))
  expect_true(all(loo$errors > 0))
})
