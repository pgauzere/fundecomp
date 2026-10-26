test_that("priority metrics use documented definitions", {
  comm <- matrix(c(2, 1, 0, 1, 0, 2), nrow = 2, byrow = TRUE,
                 dimnames = list(c("p1", "p2"), c("sp1", "sp2", "sp3")))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4))
  x <- make_fundecomp_input(comm = comm, traits = traits, time = c(1, 2))
  st <- fundecomp:::.fd_states(x)[[1L]]
  expect_equal(rarefy_community(st, "biomass", 0, "full")$value, 2 * 1 + 1 * 2)
  expect_equal(rarefy_community(st, "CWM", 0, "full")$value, 4 / 3)
  expect_equal(rarefy_community(st, "FRic", 0, "full")$value, 1)
  expect_true(all(c("FRic", "FDis", "biomass", "CWM", "variance") %in% available_fmetrics()))
})

test_that("FRic is a multivariate convex-hull volume", {
  triangle <- matrix(c(5, 2, 1), nrow = 1,
                     dimnames = list("p1", c("a", "b", "c")))
  traits_2d <- data.frame(
    species = c("a", "b", "c"),
    x = c(0, 1, 0), y = c(0, 0, 1)
  )
  x_2d <- make_fundecomp_input(comm = triangle, traits = traits_2d)
  state_2d <- fundecomp:::.fd_states(x_2d)[[1L]]
  expect_equal(
    rarefy_community(state_2d, "FRic", 0, "full", trait = c("x", "y"))$value,
    0.5, tolerance = 1e-12
  )

  tetrahedron <- matrix(c(1, 8, 2, 4), nrow = 1,
                        dimnames = list("p1", letters[1:4]))
  traits_3d <- data.frame(
    species = letters[1:4], x = c(0, 1, 0, 0),
    y = c(0, 0, 1, 0), z = c(0, 0, 0, 1)
  )
  x_3d <- make_fundecomp_input(comm = tetrahedron, traits = traits_3d)
  state_3d <- fundecomp:::.fd_states(x_3d)[[1L]]
  expect_equal(
    rarefy_community(state_3d, "FRic", 0, "full", trait = c("x", "y", "z"))$value,
    1 / 6, tolerance = 1e-12
  )

  collinear <- traits_2d
  collinear$y <- c(0, 0, 0)
  x_flat <- make_fundecomp_input(comm = triangle, traits = collinear)
  state_flat <- fundecomp:::.fd_states(x_flat)[[1L]]
  expect_equal(
    rarefy_community(state_flat, "FRic", 0, "full", trait = c("x", "y"))$value,
    0
  )
})

test_that("FRic returns zero before sampled species span trait space", {
  comm <- matrix(c(2, 1, 1), nrow = 1,
                 dimnames = list("p1", c("a", "b", "c")))
  traits <- data.frame(species = c("a", "b", "c"),
                       x = c(0, 1, 0), y = c(0, 0, 1))
  x <- make_fundecomp_input(comm = comm, traits = traits)
  state <- fundecomp:::.fd_states(x)[[1L]]
  out <- rarefy_community(state, "FRic", effort = 1:2,
                          effort_unit = "individual", trait = c("x", "y"),
                          n_rep = 20, seed = 7)
  expect_true(all(out$value >= 0, na.rm = TRUE))
  expect_false(any(is.na(out$value)))
})

test_that("individual sampling is reproducible", {
  comm <- matrix(c(3, 2, 1), nrow = 1, dimnames = list("p1", c("a", "b", "c")))
  traits <- data.frame(species = c("a", "b", "c"), mass = c(1, 2, 4))
  x <- make_fundecomp_input(comm = comm, traits = traits)
  st <- fundecomp:::.fd_states(x)[[1L]]
  a <- rarefy_community(st, "biomass", c(1, 2, 3), "individual", seed = 12)
  b <- rarefy_community(st, "biomass", c(1, 2, 3), "individual", seed = 12)
  expect_equal(a, b)
})

test_that("spatial kNN effort is available for functional curves", {
  comm <- matrix(c(2, 1, 0, 1, 0, 2, 1, 1, 1), nrow = 3, byrow = TRUE,
                 dimnames = list(paste0("p", 1:3), paste0("sp", 1:3)))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4))
  attr <- data.frame(x = c(0, 1, 0), y = c(0, 0, 1))
  x <- make_fundecomp_input(comm = comm, traits = traits, plot_attr = attr,
                    coord_names = c("x", "y"))
  st <- fundecomp:::.fd_states(x)[[1L]]
  out <- rarefy_community(st, "FDis", effort = 1:3,
                          effort_unit = "sample", spat_algo = "kNN")
  expect_true(all(is.finite(out$value)))
})
