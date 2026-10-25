sim_data <- function(sites, seed) {
  set.seed(seed)
  dat <- data.frame(time = 1:50)
  for (s in sites) {
    dat[[s]] <- rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(1) * cos(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1)
  }
  dat
}

sites  <- paste0("s", 1:10)
set.seed(1)
coords <- cbind(x = runif(10), y = runif(10))
rownames(coords) <- sites

a <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")
b <- spatfd(sim_data(sites[1:7], 2), coords, basis = "fourier", nbasis = 11, nharm = 1, name = "b")

exponential <- spatfd_model(a, model = list(a.PC1 = gstat::vgm(20, "Exp", 0.4, nugget = 1),
                                            a.PC2 = gstat::vgm(2, "Exp", 0.4)))
spherical   <- spatfd_model(a, model = list(a.PC1 = gstat::vgm(20, "Sph", 0.8, nugget = 1),
                                            a.PC2 = gstat::vgm(2, "Sph", 0.8)))
both        <- spatfd_model(a, b, model = list(a.PC1 = gstat::vgm(20, "Exp", 0.4, nugget = 1),
                                               b.PC1 = gstat::vgm(5, "Exp", 0.4, nugget = 1)))

# ================================ SALIDA ======================================
test_that("leave_one_out() returns one table per variable", {
  loo <- leave_one_out(exponential, spherical)

  expect_s3_class(loo, "SpatFDloo")
  expect_named(loo$summary, "a")
  expect_named(loo$summary$a, c("model", "method", "mean", "sd", "max"))
  expect_equal(nrow(loo$summary$a), 4)                         # 2 modelos x 2 métodos
  expect_false(is.unsorted(loo$summary$a$mean))                # Del menor al mayor error
  expect_setequal(loo$summary$a$model, c("exponential", "spherical"))

  expect_equal(nrow(loo$errors$a), 10 * 4)                     # Un error por sitio y combinación
  expect_named(loo$pred$a, c("exponential.scores", "exponential.lambda", "spherical.scores", "spherical.lambda"))
  expect_equal(dim(loo$pred$a$exponential.scores), c(50, 10))
  expect_output(print(loo), "Leave-one-out cross-validation")
})

test_that("every variable of the models is validated, each with its own models", {
  loo <- leave_one_out(exponential, both, method = "scores")
  expect_named(loo$summary, c("a", "b"))
  expect_setequal(loo$summary$a$model, c("exponential", "both"))
  expect_equal(loo$summary$b$model, "both")                    # Solo 'both' tiene a b
  expect_equal(nrow(loo$errors$b), 7)

  expect_named(leave_one_out(both, variable = "b", method = "scores")$summary, "b")
  expect_error(leave_one_out(both, variable = "c"), "Variables not in the models: c")
})

test_that("the models are labelled as in plot()", {
  loo <- leave_one_out(exponential, other = spherical, method = "scores")
  expect_setequal(loo$summary$a$model, c("exponential", "other"))
  expect_error(leave_one_out(a), "must be of class SpatFDmodel")
  expect_error(leave_one_out(exponential, method = "kriging"), "'method' must be")
})

# ============================ CADA PARTICIÓN ==================================
test_that("a site is left out of the whole analysis of its variable", {
  without <- refit_without(get_model_data(exponential)$a, "s4")

  expect_equal(rownames(without$scores), setdiff(sites, "s4"))
  expect_equal(ncol(get_pca(without)$scores), 2)
  expect_false(isTRUE(all.equal(get_pca(without)$meanfd$coefs, get_pca(a)$meanfd$coefs)))
})

test_that("the components keep the sign of the whole analysis", {
  copy    <- get_model_data(exponential)$a
  without <- refit_without(copy, "s4")
  G <- inner_product(get_pca(without)$harmonics, get_pca(copy)$harmonics, c(1, 50))
  expect_true(all(diag(G) > 0))

  # Una componente con el signo cambiado vuelve a su signo
  flipped <- without
  internal <- attr(flipped, "internal")
  internal$pca$harmonics$coefs[, 1] <- -internal$pca$harmonics$coefs[, 1]
  internal$pca$scores[, 1]          <- -internal$pca$scores[, 1]
  attr(flipped, "internal") <- internal
  flipped$scores[, "a.PC1"] <- -flipped$scores[, "a.PC1"]

  aligned <- align_signs(flipped, copy)
  expect_equal(get_pca(aligned)$scores, get_pca(without)$scores)
  expect_equal(aligned$scores, without$scores)
})

test_that("the error of a site is the prediction made without it", {
  loo <- leave_one_out(exponential, method = "scores")

  fold <- spatfd_model(refit_without(get_model_data(exponential)$a, "s4"),
                       model = get_compiled(exponential)$input)
  p    <- predict(fold, coords["s4", , drop = FALSE])
  obs  <- fda::eval.fd(1:50, get_fd(a))[, 4]

  expect_equal(unname(loo$pred$a$exponential.scores[, "s4"]), unname(p$a$pred[, 1]))
  expect_equal(loo$errors$a$error[loo$errors$a$site == "s4"], sqrt(mean((obs - p$a$pred[, 1])^2)))
})

test_that("the error is a curve in time, summarised by its root mean square", {
  fd <- fda::fd(matrix(2, 1, 1), fda::create.constant.basis(c(0, 10)))          # La curva constante 2
  e  <- error_curve(fd, matrix(c(rep(3, 6), rep(1, 5))), 0:10)
  expect_equal(e, c(rep(-1, 6), rep(1, 5)))
  expect_equal(sqrt(mean(e^2)), 1)
})

test_that("plot() draws the error curves of every model and method", {
  loo <- leave_one_out(exponential, both)
  curves <- attr(loo, "internal")$curves$a
  expect_named(curves, c("exponential.scores", "exponential.lambda", "both.scores", "both.lambda"))
  expect_equal(dim(curves$exponential.scores), c(50, 10))

  # La curva de error es la observada menos la predicha
  obs <- fda::eval.fd(1:50, get_fd(a))
  expect_equal(unname(curves$exponential.scores), unname(obs - loo$pred$a$exponential.scores))

  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(loo))
  expect_invisible(plot(loo, variable = "b"))
  expect_error(plot(loo, variable = "c"), "Variable c was not validated")
})
