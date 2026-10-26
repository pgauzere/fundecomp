test_that("mobr_path telescopes exactly at full effort", {
  comm <- matrix(c(4, 1, 0, 2, 0, 3, 1, 1, 2, 1, 2, 0), nrow = 3, byrow = TRUE,
                 dimnames = list(paste0("p", 1:3), paste0("sp", 1:4)))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4, 8))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = rep("s", 3), time = c(1, 1, 2))
  out <- get_fdelta_curves(x, list(series = "s", time = 1), list(series = "s", time = 2),
                           metric = "biomass", effort_unit = "full")
  expect_lt(max(abs(out$curves$reconstruction_error)), 1e-8)
  expect_equal(nrow(out$curves), 5)
})

test_that("Shapley sensitivity reconstructs the endpoint", {
  comm <- matrix(c(4, 1, 0, 2, 0, 3, 1, 1, 2, 1, 2, 0), nrow = 3, byrow = TRUE,
                 dimnames = list(paste0("p", 1:3), paste0("sp", 1:4)))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4, 8))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = rep("s", 3), time = c(1, 1, 2))
  out <- get_fdelta_curves(x, list(series = "s", time = 1), list(series = "s", time = 2),
                           metric = "CWM", method = "shapley", n_perm = 120)
  expect_lt(max(abs(out$curves$reconstruction_error)), 1e-8)
})

test_that("Shapley evaluates each coalition only once", {
  comm <- rbind(reference = c(a = 8, b = 5, c = 3, d = 1),
                focal = c(a = 2, b = 4, c = 6, d = 10))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4, 8))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = c("s", "s"),
                    time = 1:2)
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  counted_metric <- function(abundance, traits, trait = NULL, ...) {
    calls$n <- calls$n + 1L
    sum(abundance * traits$mass)
  }
  out <- get_fdelta_curves(
    x, list(series = "s", time = 1), list(series = "s", time = 2),
    metric = counted_metric, method = "shapley", n_perm = 120
  )
  expect_lt(max(abs(out$curves$reconstruction_error)), 1e-8)
  expect_lte(calls$n, 34L) # 32 coalitions plus two explicit endpoint checks
})

test_that("nested state construction enforces N and SAD margins", {
  comm <- matrix(c(8, 2, 0, 1, 3, 6), nrow = 2, byrow = TRUE,
                 dimnames = list(c("p1", "p2"), c("sp1", "sp2", "sp3")))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = c("s", "s"), time = 1:2)
  pair <- fundecomp:::.fd_prepare_pair(list(series = "s", time = 1),
                                     list(series = "s", time = 2), x)
  state <- fundecomp:::.fd_subset_state(pair, c("N", "SAD"), x)
  expect_equal(sum(state$comm), sum(pair$focal$comm), tolerance = 1e-8)
  expect_equal(colSums(state$comm) / sum(state$comm),
               colSums(pair$focal$comm) / sum(pair$focal$comm), tolerance = 1e-8)
  expect_equal(rowSums(state$comm) / sum(state$comm),
               rowSums(pair$ref$comm) / sum(pair$ref$comm), tolerance = 1e-8)
})

test_that("sample effort remains valid when dates have different site counts", {
  comm <- matrix(c(4, 1, 0, 2,
                   0, 3, 1, 1,
                   2, 0, 2, 1,
                   1, 2, 0, 3,
                   3, 1, 1, 0), nrow = 5, byrow = TRUE,
                 dimnames = list(paste0("p", 1:5), paste0("sp", 1:4)))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4, 8))
  x <- make_fundecomp_input(
    comm = comm, traits = traits, series = rep("s", 5),
    time = c(1, 1, 2, 2, 2),
    sample = paste0("p", 1:5),
    plot_attr = data.frame(longitude = c(0, 1, 0, 1, 2),
                           latitude = c(0, 0, 1, 1, 0)),
    coord_names = c("longitude", "latitude"), latlong = TRUE
  )
  out <- get_fdelta_curves(
    x, list(series = "s", time = 1), list(series = "s", time = 2),
    metric = "FDis", effort = c(1, 2), effort_unit = "sample",
    spat_algo = "kNN"
  )
  expect_true(all(is.finite(out$curves$estimate)))
  expect_lt(max(abs(out$curves$reconstruction_error)), 1e-8)
})
