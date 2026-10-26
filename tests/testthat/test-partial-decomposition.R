partial_example <- function() {
  comm <- rbind(r1 = c(a = 5, b = 3, c = 1, d = 0),
                r2 = c(a = 2, b = 1, c = 3, d = 2),
                f1 = c(a = 1, b = 2, c = 2, d = 5),
                f2 = c(a = 0, b = 1, c = 4, d = 3))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4, 8),
                       pace = c(5, 1, 4, 2))
  make_fundecomp_input(comm = comm, traits = traits,
    series = rep("s", 4), time = c(1, 1, 2, 2))
}

partial_run <- function(x, components, method = "mobr_path", metric = "CWM",
                        trait = "mass", ta_axis = NULL,
                        effort_unit = "full", n_rep = 1, seed = 27,
                        tie_method = "average", ...) {
  decompose_fmetric_change(x, list(series = "s", time = 1),
    list(series = "s", time = 2), metric = metric, trait = trait, ta_axis = ta_axis,
    effort_unit = effort_unit, n_rep = n_rep, seed = seed, tie_method = tie_method,
    components = components, method = method, ...)
}

test_that("selected components reconstruct their conditional contrast", {
  x <- partial_example()
  for (selected in list("N", c("N", "SAD"), c("STD", "TA"),
                        c("N", "aggregation", "SAD"))) {
    result <- partial_run(x, selected)
    expect_equal(sum(result$curves$estimate), result$curves$partial_change[1],
                 tolerance = 1e-12)
    expect_equal(result$curves$partial_change[1],
                 result$curves$partial_counterfactual_metric[1] -
                   result$curves$reference_metric[1], tolerance = 1e-12)
    expect_equal(result$curves$observed_total_change[1],
                 result$curves$observed_focal_metric[1] -
                   result$curves$reference_metric[1], tolerance = 1e-12)
    expect_equal(result$curves$reconstruction_residual[1],
                 result$curves$partial_change[1] - sum(result$curves$estimate),
                 tolerance = 1e-12)
  }
  partial <- partial_run(x, c("N", "SAD"))
  expect_equal(partial$components, c("N", "SAD"))
  expect_equal(partial$curves$remaining_difference[1],
               partial$curves$observed_total_change[1] - partial$curves$partial_change[1],
               tolerance = 1e-12)
  expect_true(validate_decomposition(partial)$passed)
})

test_that("selecting every component reproduces the default full decomposition", {
  x <- partial_example()
  default <- partial_run(x, NULL)
  selected <- partial_run(x, fundecomp:::.fd_components_for(x))
  expect_equal(selected$curves$estimate, default$curves$estimate, tolerance = 0)
  expect_equal(selected$curves$reconstruction_error,
               default$curves$reconstruction_error, tolerance = 0)
  expect_equal(selected$curves$partial_change,
               selected$curves$observed_total_change, tolerance = 0)
})

test_that("restricted Shapley averages the selected-component game", {
  x <- partial_example()
  result <- partial_run(x, c("N", "SAD"), method = "shapley", n_rep = 2)
  expect_equal(length(result$permutations), 2L)
  expect_equal(sum(result$curves$estimate), result$curves$partial_change[1],
               tolerance = 1e-12)
  expect_true(all(result$curves$selected_components[[1]] == c("N", "SAD")))
  exact <- shapley_decomposition(partial_run(x, c("N", "SAD")), exact = TRUE)
  expect_equal(length(exact$permutations), 2L)
  expect_equal(sum(exact$curves$estimate), exact$curves$partial_change[1],
               tolerance = 1e-12)
})

test_that("partial decomposition can select ITV while holding other components fixed", {
  comm <- rbind(reference = c(a = 10L, b = 10L), focal = c(a = 10L, b = 10L))
  traits <- data.frame(species = c("a", "b"), mass = c(0, 2))
  observations <- data.frame(species = rep(c("a", "b"), each = 4L),
    series = "s", time = rep(c(1, 1, 2, 2), 2L),
    mass = c(-1, 1, -2, 2, 1, 3, 0, 4))
  x <- make_fundecomp_input(comm = comm, traits = traits,
    series = rep("s", 2L), time = 1:2, trait_observations = observations)
  result <- decompose_fmetric_change(
    x, list(series = "s", time = 1), list(series = "s", time = 2),
    metric = "variance", components = "ITV", effort_unit = "full",
    n_rep = 1, seed = 5
  )
  expect_equal(result$curves$estimate, 3, tolerance = 1e-12)
  expect_equal(result$curves$partial_change, 3, tolerance = 1e-12)
  expect_equal(result$curves$observed_total_change, 3, tolerance = 1e-12)
})

test_that("partial decompositions work in temporal contrasts and multivariate metrics", {
  x <- partial_example()
  temporal <- decompose_fmetric_ts(x, metric = "FDis", components = c("STD", "TA"),
    trait = c("mass", "pace"), ta_axis = "mass", effort_unit = "full",
    n_rep = 1, seed = 9, tie_method = "average")
  expect_equal(unique(temporal$curves$component), c("traits", "TA"))
  expect_lt(max(abs(temporal$curves$reconstruction_residual)), 1e-12)
  expect_equal(temporal$curves$estimate,
    partial_run(x, c("STD", "TA"), metric = "FDis", trait = c("mass", "pace"),
                ta_axis = "mass")$curves$estimate,
    tolerance = 1e-12)
})

test_that("partial results repeat exactly under a fixed seed", {
  x <- partial_example()
  a <- partial_run(x, c("N", "SAD"), effort = 1:2, effort_unit = "sample")
  b <- partial_run(x, c("N", "SAD"), effort = 1:2, effort_unit = "sample")
  expect_equal(a$curves, b$curves, tolerance = 0)
})
