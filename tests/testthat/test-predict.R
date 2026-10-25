sim_data <- function(sites, seed, time = 1:50) {
  set.seed(seed)
  dat <- data.frame(time = time)
  for (s in sites) {
    dat[[s]] <- rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(1) * cos(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1)
  }
  dat
}

sites  <- paste0("s", 1:12)
set.seed(1)
coords <- cbind(x = runif(12), y = runif(12))
rownames(coords) <- sites

a <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")
b <- spatfd(sim_data(sites[1:8], 2), coords, basis = "fourier", nbasis = 11, nharm = 1, name = "b")

direct <- list(a.PC1 = gstat::vgm(20, "Exp", 0.4, nugget = 1),
               a.PC2 = gstat::vgm(2, "Sph", 0.6))
m <- spatfd_model(a, model = direct)

new <- rbind(cbind(x = c(0.2, 0.7, 0.5), y = c(0.5, 0.5, 0.1)))
rownames(new) <- c("n1", "n2", "n3")
new_sp <- as.data.frame(new)
sp::coordinates(new_sp) <- ~ x + y

# ============================ CONTRA gstat ====================================
test_that("the scores are kriged as gstat does", {
  k <- krige_scores(m, new)
  g <- to_gstat(a)

  for (id in c("a.PC1", "a.PC2")) {
    ref <- gstat::krige(stats::as.formula(paste(id, "~ 1")), g$data[[id]]$data, new_sp,
                        model = direct[[id]], debug.level = 0)
    expect_equal(unname(k$pred[, id]), ref$var1.pred, info = id)
    expect_equal(unname(k$variance[, id]), ref$var1.var, info = id)
  }
})

test_that("the scores are cokriged as gstat does, with sites of their own", {
  g   <- to_gstat(a, b, pcs = c(a = 1, b = 1))
  fit <- suppressWarnings(gstat::fit.lmc(gstat::variogram(g), g, gstat::vgm(1, "Exp", 0.5, nugget = 0.1),
                                         fit.method = 6, correct.diagonal = 1.01))
  k   <- krige_scores(spatfd_model(a, b, model = fit), new)
  ref <- suppressMessages(predict(fit, new_sp, debug.level = 0))

  for (id in c("a.PC1", "b.PC1")) {
    expect_equal(unname(k$pred[, id]), ref[[paste0(id, ".pred")]], info = id)
    expect_equal(unname(k$variance[, id]), ref[[paste0(id, ".var")]], info = id)
  }
})

# ============================== CURVAS ========================================
test_that("predict() rebuilds the curves from the predicted scores", {
  p   <- predict(m, new)
  k   <- krige_scores(m, new)
  pca <- get_pca(get_model_data(m)$a)

  expect_s3_class(p, "SpatFDpred")
  expect_named(p, "a")
  expect_equal(dim(p$a$pred), c(50, 3))
  expect_equal(colnames(p$a$pred), c("n1", "n2", "n3"))

  curve <- fda::eval.fd(1:50, pca$meanfd) + fda::eval.fd(1:50, pca$harmonics) %*% k$pred["n2", ]
  expect_equal(unname(p$a$pred[, "n2"]), as.vector(curve))
  expect_equal(unname(p$a$var), unname(rowSums(k$variance)))
})

test_that("without a nugget an observed site is reproduced", {
  exact <- spatfd_model(a, model = list(a.PC1 = gstat::vgm(20, "Exp", 0.4), a.PC2 = gstat::vgm(2, "Sph", 0.6)))
  site  <- coords["s3", , drop = FALSE]

  k <- krige_scores(exact, site)
  expect_equal(unname(k$pred[1, ]), unname(get_pca(a)$scores[3, ]))
  expect_equal(unname(k$variance[1, ]), c(0, 0), tolerance = 1e-8)

  # El método lambda devuelve la curva observada completa
  p <- predict(exact, site, method = "lambda")
  expect_equal(unname(p$a$pred[, 1]), as.vector(fda::eval.fd(1:50, get_fd(a))[, 3]))
  expect_equal(unname(p$a$var), 0, tolerance = 1e-8)
})

test_that("the method lambda predicts every variable of a cokriging model", {
  g   <- to_gstat(a, b, pcs = c(a = 1, b = 1))
  fit <- suppressWarnings(gstat::fit.lmc(gstat::variogram(g), g, gstat::vgm(1, "Exp", 0.5, nugget = 0.1),
                                         fit.method = 6, correct.diagonal = 1.01))
  p <- predict(spatfd_model(a, b, model = fit), new, method = "lambda")

  expect_named(p, c("a", "b"))
  expect_true(all(p$a$var > 0))
  expect_equal(dim(p$b$pred), c(50, 3))
})

test_that("the method lambda needs times in common", {
  late <- spatfd(sim_data(sites, 3, time = 101:150), coords, basis = "fourier", nbasis = 11, nharm = 1,
                 name = "late")
  both <- spatfd_model(a, late, model = list(a.PC1 = gstat::vgm(1, "Exp", 0.4), late.PC1 = gstat::vgm(1, "Exp", 0.4)))

  expect_no_error(predict(both, new))
  expect_error(predict(both, new, method = "lambda"), "no time in common")
})

# ============================== ARGUMENTOS ====================================
test_that("variable selects what is returned", {
  g   <- to_gstat(a, b, pcs = c(a = 1, b = 1))
  fit <- suppressWarnings(gstat::fit.lmc(gstat::variogram(g), g, gstat::vgm(1, "Exp", 0.5, nugget = 0.1),
                                         fit.method = 6, correct.diagonal = 1.01))
  ab  <- spatfd_model(a, b, model = fit)

  expect_named(predict(ab, new), c("a", "b"))
  expect_named(predict(ab, new, variable = "b"), "b")
  expect_equal(predict(ab, new, variable = "b")$b, predict(ab, new)$b)
  expect_error(predict(ab, new, variable = "c"), "Variables not in the model: c")
})

test_that("newcoords are matched by name or by position", {
  unnamed <- unname(new)
  p <- predict(m, unnamed)
  expect_equal(colnames(p$a$pred), c("new1", "new2", "new3"))
  expect_equal(unname(p$a$pred), unname(predict(m, new)$a$pred))
  expect_equal(predict(m, new[, c("y", "x")])$a$pred, predict(m, new)$a$pred)
  expect_error(predict(m, cbind(lon = 1, lat = 2)), "must have the coordinates of the model")
})

test_that("the rows of the prediction are the times of the data", {
  daily <- spatfd(sim_data(sites, 1, time = as.Date("2026-01-01") + 0:49), coords, basis = "fourier",
                  nbasis = 11, nharm = 2, name = "a")
  p <- predict(spatfd_model(daily, model = direct), new)
  expect_equal(rownames(p$a$pred)[1:2], c("2026-01-01", "2026-01-02"))
  expect_equal(rownames(predict(m, new)$a$pred)[1:2], c("1", "2"))
})

# ================================= MAPAS ======================================
test_that("predict() works over a map", {
  skip_if_not_installed("sf")
  square <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))))
  p <- predict(m, square, grid_size = 100)
  expect_equal(ncol(p$a$pred), nrow(attr(p, "newcoords")))
  expect_false(is.null(attr(p, "map")))

  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(p, time = 10))
  expect_invisible(plot(p, mode = "curves"))
})

# ================================ MÉTODOS =====================================
test_that("print() and plot() work", {
  p <- predict(m, new)
  expect_output(print(p), "Predicted curves at 3 sites")

  daily <- spatfd(sim_data(sites, 1, time = as.Date("2026-01-01") + 0:49), coords, basis = "fourier",
                  nbasis = 11, nharm = 2, name = "a")
  pd <- predict(spatfd_model(daily, model = direct), new)

  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(p))
  expect_invisible(plot(pd, col = "red"))
  expect_error(plot(p, mode = "map"), "only available for predictions made over a map")
})
