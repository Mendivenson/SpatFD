sim_data <- function(sites, seed) {
  set.seed(seed)
  dat <- data.frame(time = 1:50)
  for (s in sites) {
    dat[[s]] <- rnorm(1, 5) * sin(2 * pi * (1:50) / 50) + rnorm(1) * cos(2 * pi * (1:50) / 50) + rnorm(50, sd = 0.1)
  }
  dat
}

# Estaciones solo en la mitad izquierda del cuadrado unitario
sites  <- paste0("s", 1:8)
set.seed(1)
coords <- cbind(x = runif(8, 0, 0.4), y = runif(8))
rownames(coords) <- sites

a <- spatfd(sim_data(sites, 1), coords, basis = "fourier", nbasis = 11, nharm = 2, name = "a")
model <- spatfd_model(a, model = list(a.PC1 = gstat::vgm(20, "Exp", 0.3, nugget = 1),
                                      a.PC2 = gstat::vgm(2, "Exp", 0.3)))

targets    <- as.matrix(expand.grid(x = seq(0, 1, 0.1), y = seq(0, 1, 0.1)))
candidates <- as.matrix(expand.grid(x = seq(0.05, 0.95, 0.15), y = seq(0.05, 0.95, 0.15)))

# ============================= LO QUE SE MINIMIZA =============================
test_that("the total variance is the sum of the prediction variances at the targets", {
  for (method in c("scores", "lambda")) {
    systems <- get_compiled(model)[[method]]
    if (method == "lambda") systems <- systems$systems

    expected <- sum(predict(model, targets, method = method)$a$var)
    expect_equal(design_variance(systems, candidates[0, , drop = FALSE], check_newcoords(targets, c("x", "y"))),
                 expected, info = method)
  }
})

test_that("optimal_design() reports the variance of the design it returns", {
  design <- optimal_design(model, k = 2, targets = targets, candidates = candidates)

  expect_s3_class(design, "SpatFDdesign")
  expect_equal(dim(design$new_stations), c(2, 2))
  expect_equal(design$variance$step, c("existing stations", "+ new1", "+ new2", "after exchanges"))
  expect_equal(design$variance$mean, design$variance$total / nrow(targets))

  # La varianza baja con cada estación y los intercambios no la suben
  expect_false(is.unsorted(rev(design$variance$total)))

  # Es la varianza de predecir con las estaciones nuevas, en el sistema del modelo
  systems <- get_compiled(model)$scores
  expect_equal(design_variance(systems, design$new_stations, check_newcoords(targets, c("x", "y"))),
               design$variance$total[4])
})

test_that("the fast variance of a design is the exact one", {
  t <- check_newcoords(targets, c("x", "y"))
  for (method in c("scores", "lambda")) {
    systems <- get_compiled(model)[[method]]
    if (method == "lambda") systems <- systems$systems
    cache <- design_cache(systems, candidates, t)

    for (design in list(integer(0), 5L, c(3L, 40L), c(1L, 20L, 33L, 49L))) {
      exact <- design_variance(systems, candidates[design, , drop = FALSE], t)
      expect_equal(cache$value(design), exact, info = paste(method, toString(design)))
    }

    # Agregar cada candidato a la vez da lo mismo que evaluarlos uno por uno
    for (design in list(integer(0), c(3L, 40L))) {
      added <- cache$add(design)
      for (i in c(1L, 17L, 49L)) {
        exact <- design_variance(systems, candidates[c(design, i), , drop = FALSE], t)
        expect_equal(added[i], exact, info = paste(method, toString(design), i))
      }
      expect_true(all(is.infinite(added[design])))
    }
  }
})

test_that("candidates in groups of several fields are added one by one", {
  # Un variograma cruzado entre las dos componentes deja un grupo de dos campos
  cross <- spatfd_model(a, model = list(a.PC1 = gstat::vgm(20, "Exp", 0.3, nugget = 1),
                                        a.PC2 = gstat::vgm(2, "Exp", 0.3),
                                        a.PC1.a.PC2 = gstat::vgm(1, "Exp", 0.3)))
  t       <- check_newcoords(targets, c("x", "y"))
  systems <- get_compiled(cross)$scores
  cache   <- design_cache(systems, candidates, t)

  expect_equal(cache$value(c(3L, 40L)), design_variance(systems, candidates[c(3, 40), ], t))
  expect_equal(cache$add(3L)[40], design_variance(systems, candidates[c(3, 40), ], t))
})

# ============================== LA BÚSQUEDA ===================================
test_that("the first station goes where the prediction is worst", {
  near <- cbind(x = 0.2, y = 0.5)
  far  <- cbind(x = 0.9, y = 0.5)
  design <- optimal_design(model, k = 1, targets = targets, candidates = rbind(near, far))
  expect_equal(unname(design$new_stations[1, ]), c(0.9, 0.5))
})

test_that("the design is deterministic and never repeats a station", {
  d1 <- optimal_design(model, k = 3, targets = targets, candidates = candidates)
  d2 <- optimal_design(model, k = 3, targets = targets, candidates = candidates)
  expect_equal(d1, d2)
  expect_false(any(duplicated(d1$new_stations)))
})

test_that("random candidates are inside the area and can be repeated", {
  set.seed(5)
  d1 <- optimal_design(model, k = 1, targets = targets, candidate_size = 30)
  set.seed(5)
  d2 <- optimal_design(model, k = 1, targets = targets, candidate_size = 30)
  expect_equal(d1$new_stations, d2$new_stations)

  candidates <- attr(d1, "candidates")
  expect_equal(nrow(candidates), 30)
  expect_true(all(candidates >= 0 & candidates <= 1))
})

test_that("targets can be a map", {
  skip_if_not_installed("sf")
  square <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))))
  set.seed(2)
  design <- optimal_design(model, k = 1, targets = square, target_size = 50, candidate_size = 20)
  expect_false(is.null(attr(design, "map")))
  expect_true(all(attr(design, "candidates") >= 0 & attr(design, "candidates") <= 1))

  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(design))
})

# ============================== ARGUMENTOS ====================================
test_that("optimal_design() validates its arguments", {
  b  <- spatfd(sim_data(sites, 2), coords, basis = "fourier", nbasis = 11, nharm = 1, name = "b")
  ab <- spatfd_model(a, b, model = list(a.PC1 = gstat::vgm(1, "Exp", 0.3), b.PC1 = gstat::vgm(1, "Exp", 0.3)))

  expect_error(optimal_design(ab, k = 1, targets = targets), "single variable")
  expect_error(optimal_design(model, k = 0, targets = targets), "'k' must be")
  expect_error(optimal_design(model, k = 3, targets = targets, candidates = candidates[1:2, ]),
               "cannot be greater than the number of candidates \\(2\\)")

  # Un candidato sobre una estación existente se descarta
  design <- optimal_design(model, k = 1, targets = targets, candidates = rbind(coords["s1", ], c(0.9, 0.5)))
  expect_equal(nrow(attr(design, "candidates")), 1)
})

test_that("progress bars show the candidates tried in each stage", {
  expect_message(
    expect_output(optimal_design(model, k = 2, targets = targets, candidates = candidates, progress = TRUE),
                  "100%"),
    "Adding 2 stations: candidates tried")
  expect_silent(optimal_design(model, k = 2, targets = targets, candidates = candidates, progress = FALSE))
  expect_error(optimal_design(model, k = 2, targets = targets, progress = "yes"), "'progress' must be")
})

test_that("print() and plot() work", {
  design <- optimal_design(model, k = 2, targets = targets, candidates = candidates)
  expect_output(print(design), "Optimal design")

  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(design))
})
