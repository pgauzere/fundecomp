test_that("CWV is the abundance-weighted population variance", {
  comm <- matrix(c(2, 1), nrow = 1, dimnames = list("p1", c("a", "b")))
  traits <- data.frame(species = c("a", "b"), mass = c(1, 3))
  x <- make_fundecomp_input(comm = comm, traits = traits)
  expect_equal(rarefy_community(x, metric = "CWV", effort = 0,
                                effort_unit = "full", trait = "mass")$value,
               8 / 9, tolerance = 1e-12)
})

test_that("single-trait metrics require an explicit trait in multivariate input", {
  comm <- matrix(c(1, 2), nrow = 1, dimnames = list("p1", c("a", "b")))
  traits <- data.frame(species = c("a", "b"), mass = c(1, 2), pace = c(3, 4))
  x <- make_fundecomp_input(comm = comm, traits = traits)
  expect_error(rarefy_community(x, metric = "CWM", effort = 0,
                                effort_unit = "full"), "supply `trait`")
  expect_error(rarefy_community(x, metric = "CWV", effort = 0,
                                effort_unit = "full"), "supply `trait`")
  expect_error(rarefy_community(x, metric = "trait_sum", effort = 0,
                                effort_unit = "full"), "supply `trait`")
})

test_that("generic state labels support endpoint comparison", {
  comm <- rbind(control = c(a = 2, b = 1), treatment = c(a = 1, b = 2))
  traits <- data.frame(species = c("a", "b"), mass = c(1, 4))
  x <- make_fundecomp_input(comm = comm, traits = traits, state = c("control", "treatment"))
  out <- decompose_fmetric_change(
    x, reference = list(series = "series1", state = "control"),
    focal = list(series = "series1", state = "treatment"),
    metric = "CWM", trait = "mass", effort_unit = "full", n_rep = 1, seed = 3
  )
  expect_lt(max(abs(out$curves$reconstruction_error)), 1e-12)
  expect_false("pool" %in% names(out))
})

test_that("rarefaction refuses implicit selection from a multi-state object", {
  comm <- rbind(t1 = c(a = 2, b = 1), t2 = c(a = 1, b = 2))
  traits <- data.frame(species = c("a", "b"), mass = c(1, 4))
  x <- make_fundecomp_input(comm = comm, traits = traits, state = c("t1", "t2"))
  expect_error(rarefy_community(x, metric = "CWM", effort = 0,
                                effort_unit = "full", trait = "mass"),
               "must be specified")
})

test_that("Shapley is a secondary exact analysis of a path object", {
  comm <- rbind(r = c(a = 4, b = 2, c = 1), f = c(a = 1, b = 2, c = 4))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = c("s", "s"), time = 1:2)
  path <- decompose_fmetric_change(
    x, list(series = "s", state = 1), list(series = "s", state = 2),
    metric = "CWM", trait = "mass", effort_unit = "full", n_rep = 1, seed = 8
  )
  shapley <- shapley_decomposition(path)
  expect_equal(length(shapley$permutations), 120L)
  expect_true(all(shapley$curves$interval_type == "structural_path_sensitivity"))
  expect_lt(max(abs(shapley$curves$reconstruction_error)), 1e-12)
})
