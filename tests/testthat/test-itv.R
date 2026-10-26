itv_example <- function() {
  comm <- rbind(reference = c(a = 10L, b = 10L),
                focal = c(a = 10L, b = 10L))
  traits <- data.frame(species = c("a", "b"), mass = c(0, 2))
  observations <- data.frame(
    species = rep(c("a", "b"), each = 4L),
    series = "s",
    time = rep(c(1, 1, 2, 2), 2L),
    mass = c(-1, 1, -2, 2, 1, 3, 0, 4)
  )
  make_fundecomp_input(
    comm = comm, traits = traits, series = rep("s", 2L), time = 1:2,
    trait_observations = observations
  )
}

test_that("date-specific trait observations create an ITV component", {
  x <- itv_example()
  expect_true(x$has_itv)
  expect_equal(fundecomp:::.fd_components_for(x),
               c("N", "aggregation", "SAD", "traits", "ITV", "TA"))
  states <- fundecomp:::.fd_states(x)
  expect_equal(states[[1L]]$traits["a", "mass"], 0)
  expect_equal(states[[2L]]$traits["a", "mass"], 0)
  expect_equal(as.numeric(states[[1L]]$itv$a[, "mass"]), c(-1, 1))
  expect_equal(as.numeric(states[[2L]]$itv$a[, "mass"]), c(-2, 2))
  expect_equal(
    rarefy_community(states[[1L]], "FRic", 0, "full", trait = "mass")$value,
    4, tolerance = 1e-12
  )
})

test_that("an isolated within-species trait change is assigned to ITV", {
  x <- itv_example()
  reference <- list(series = "s", time = 1)
  focal <- list(series = "s", time = 2)
  path <- get_fdelta_curves(
    x, reference, focal, metric = "variance", effort_unit = "full",
    method = "mobr_path", seed = 4
  )
  estimates <- stats::setNames(path$curves$estimate, path$curves$component)
  expect_equal(unname(estimates["ITV"]), 3, tolerance = 1e-12)
  expect_lt(max(abs(estimates[setdiff(names(estimates), "ITV")])), 1e-12)
  expect_lt(max(abs(path$curves$reconstruction_error)), 1e-12)

  shapley <- get_fdelta_curves(
    x, reference, focal, metric = "variance", effort_unit = "full",
    method = "shapley", n_perm = 720L, seed = 4
  )
  estimates <- stats::setNames(shapley$curves$estimate, shapley$curves$component)
  expect_equal(unname(estimates["ITV"]), 3, tolerance = 1e-12)
  expect_lt(max(abs(estimates[setdiff(names(estimates), "ITV")])), 1e-12)
  expect_lt(max(abs(shapley$curves$reconstruction_error)), 1e-12)
  expect_equal(length(shapley$permutations), 720L)
})

test_that("scoped trait observations require coverage or an explicit fallback", {
  comm <- rbind(reference = c(a = 2L, b = 2L), focal = c(a = 2L, b = 2L))
  traits <- data.frame(species = c("a", "b"), mass = c(0, 2))
  observations <- data.frame(
    species = c("a", "b", "a"), series = "s", time = c(1, 1, 2),
    mass = c(0, 2, 0)
  )
  x <- make_fundecomp_input(comm = comm, traits = traits, series = rep("s", 2L),
                    time = 1:2, trait_observations = observations)
  expect_error(fundecomp:::.fd_states(x), "No trait observations")

  fallback <- make_fundecomp_input(
    comm = comm, traits = traits, series = rep("s", 2L), time = 1:2,
    trait_observations = observations, itv_na = "species_mean"
  )
  states <- fundecomp:::.fd_states(fallback)
  expect_equal(unname(states[[2L]]$itv$b), matrix(0, nrow = 1L, ncol = 1L))
})

test_that("inputs without trait observations retain the five-component API", {
  comm <- rbind(reference = c(a = 2L, b = 1L), focal = c(a = 1L, b = 2L))
  traits <- data.frame(species = c("a", "b"), mass = c(0, 2))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = rep("s", 2L), time = 1:2)
  expect_false(x$has_itv)
  expect_equal(fundecomp:::.fd_components_for(x),
               c("N", "aggregation", "SAD", "traits", "TA"))
})

test_that("time-series decomposition propagates the ITV component", {
  out <- decompose_fmetric_ts(
    itv_example(), metric = "variance", time_reference = "first",
    effort_unit = "full", seed = 11
  )
  expect_true("ITV" %in% out$curves$component)
  itv_rows <- out$curves[out$curves$component == "ITV", , drop = FALSE]
  expect_equal(itv_rows$estimate, 3, tolerance = 1e-12)
  expect_lt(max(abs(itv_rows$reconstruction_error)), 1e-12)
})

test_that("multivariate FRic can isolate a change in ITV", {
  comm <- rbind(reference = c(a = 10L, b = 10L),
                focal = c(a = 10L, b = 10L))
  traits <- data.frame(species = c("a", "b"), x = c(0, 2), y = c(0, 0))
  observations <- data.frame(
    species = rep(c("a", "b"), each = 4L), series = "s",
    time = rep(c(1, 1, 2, 2), 2L),
    x = c(-1, 1, -2, 2, 1, 3, 0, 4),
    y = c(0, 0, -1, 1, -1, 1, -2, 2)
  )
  x <- make_fundecomp_input(
    comm = comm, traits = traits, series = rep("s", 2L), time = 1:2,
    trait_observations = observations
  )
  out <- get_fdelta_curves(
    x, list(series = "s", time = 1), list(series = "s", time = 2),
    metric = "FRic", trait = c("x", "y"), ta_axis = "x", effort_unit = "full"
  )
  estimates <- stats::setNames(out$curves$estimate, out$curves$component)
  expect_equal(unname(estimates["ITV"]), 3, tolerance = 1e-12)
  expect_lt(max(abs(estimates[setdiff(names(estimates), "ITV")])), 1e-12)
})
