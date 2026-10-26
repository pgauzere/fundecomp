tie_example <- function(reference, focal, mass, sample = NULL, time = c(1, 2)) {
  comm <- rbind(reference = reference, focal = focal)
  species <- paste0("sp", seq_along(reference))
  colnames(comm) <- species
  traits <- data.frame(species = species, mass = mass)
  make_fundecomp_input(
    comm = comm, traits = traits, series = rep("s", nrow(comm)),
    time = time, sample = sample
  )
}

renamed_tie_example <- function(x, permutation) {
  comm <- x$comm[, permutation, drop = FALSE]
  new_species <- paste0("renamed_", seq_along(permutation))
  colnames(comm) <- new_species
  traits <- data.frame(
    species = new_species,
    mass = x$traits$mass[permutation]
  )
  make_fundecomp_input(
    comm = comm, traits = traits, series = x$meta$series,
    time = x$meta$time, sample = x$meta$sample
  )
}

test_that("tie-averaged path is invariant to species labels", {
  x <- tie_example(
    reference = c(10, 10, 5, 5, 1, 1),
    focal = c(5, 5, 10, 10, 1, 1),
    mass = c(1, 10, 100, 1000, 10000, 100000)
  )
  x_renamed <- renamed_tie_example(x, c(4, 1, 6, 2, 5, 3))
  args <- list(
    reference = list(series = "s", time = 1),
    focal = list(series = "s", time = 2),
    metric = "biomass", effort_unit = "full",
    tie_method = "average", seed = 91
  )
  out <- do.call(get_fdelta_curves, c(list(x = x), args))
  out_renamed <- do.call(get_fdelta_curves, c(list(x = x_renamed), args))

  expect_equal(out$curves$estimate, out_renamed$curves$estimate, tolerance = 1e-12)
  expect_equal(out$curves$lower, out_renamed$curves$lower, tolerance = 1e-12)
  expect_equal(out$curves$upper, out_renamed$curves$upper, tolerance = 1e-12)
  expect_equal(out$curves$reconstruction_error,
               out_renamed$curves$reconstruction_error, tolerance = 1e-12)
  expect_true(out$diagnostics$tie$states[[2]]$tied)
  expect_true(out$diagnostics$tie$states[[2]]$exact)
})

test_that("strictly ordered abundance vectors retain previous estimates", {
  x <- tie_example(
    reference = c(8, 4, 2, 1),
    focal = c(7, 5, 2, 1),
    mass = c(1, 2, 4, 8)
  )
  args <- list(
    reference = list(series = "s", time = 1),
    focal = list(series = "s", time = 2),
    metric = "biomass", effort_unit = "full", seed = 17
  )
  averaged <- do.call(get_fdelta_curves,
                      c(list(x = x, tie_method = "average"), args))
  legacy <- NULL
  expect_warning(
    legacy <- do.call(get_fdelta_curves,
                      c(list(x = x, tie_method = "legacy"), args)),
    "label-dependent"
  )
  expect_equal(averaged$curves$estimate, legacy$curves$estimate, tolerance = 1e-12)
  expect_equal(averaged$curves$lower, legacy$curves$lower, tolerance = 1e-12)
  expect_true(all(vapply(averaged$diagnostics$tie$states,
                         function(z) !isTRUE(z$tied), logical(1))))

  shapley_average <- do.call(get_fdelta_curves,
                             c(list(x = x, tie_method = "average",
                                    method = "shapley", n_perm = 120L), args))
  shapley_legacy <- NULL
  expect_warning(
    shapley_legacy <- do.call(
      get_fdelta_curves,
      c(list(x = x, tie_method = "legacy",
             method = "shapley", n_perm = 120L), args)
    ),
    "label-dependent"
  )
  expect_equal(shapley_average$curves$estimate,
               shapley_legacy$curves$estimate, tolerance = 1e-12)
})

test_that("a pure T-A change with ties does not leak into SAD or STD", {
  x <- tie_example(
    reference = c(10, 10, 5, 5),
    focal = c(5, 5, 10, 10),
    mass = c(1, 2, 100, 101)
  )
  out <- get_fdelta_curves(
    x, list(series = "s", time = 1), list(series = "s", time = 2),
    metric = "biomass", effort_unit = "full", tie_method = "average", seed = 2
  )
  estimates <- stats::setNames(out$curves$estimate, out$curves$component)
  expect_equal(unname(estimates[c("N", "aggregation", "SAD", "traits")]),
               rep(0, 4), tolerance = 1e-12)
  expect_equal(unname(estimates["TA"]), 990, tolerance = 1e-12)
  expect_lt(max(abs(out$curves$reconstruction_error)), 1e-12)
})

test_that("complete abundance ties have a label-invariant, averaged estimand", {
  x <- tie_example(
    reference = rep(5, 4), focal = rep(5, 4),
    mass = c(1, 10, 100, 1000)
  )
  x_renamed <- renamed_tie_example(x, c(3, 1, 4, 2))
  out <- get_fdelta_curves(
    x, list(series = "s", time = 1), list(series = "s", time = 2),
    metric = "CWM", effort_unit = "full", tie_method = "average", seed = 11
  )
  out_renamed <- get_fdelta_curves(
    x_renamed, list(series = "s", time = 1), list(series = "s", time = 2),
    metric = "CWM", effort_unit = "full", tie_method = "average", seed = 11
  )
  expect_equal(out$curves$estimate, out_renamed$curves$estimate, tolerance = 1e-12)
  expect_equal(out$curves$estimate, rep(0, 5), tolerance = 1e-12)
  expect_equal(out$diagnostics$tie$states[[2]]$total, factorial(4)^2)
  expect_true(out$diagnostics$tie$states[[2]]$exact)
})

test_that("partial sample effort remains invariant with spatially distinct tied species", {
  comm <- rbind(
    r1 = c(5, 5, 0, 0), r2 = c(0, 0, 3, 3),
    f1 = c(3, 3, 0, 0), f2 = c(0, 0, 5, 5)
  )
  colnames(comm) <- paste0("sp", 1:4)
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 100, 101))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = rep("s", 4),
                    time = c(1, 1, 2, 2), sample = rownames(comm))
  x_renamed <- renamed_tie_example(x, c(3, 1, 4, 2))
  args <- list(
    reference = list(series = "s", time = 1),
    focal = list(series = "s", time = 2), metric = "CWM",
    effort = 1:2, effort_unit = "sample", tie_method = "average",
    n_perm = 20, seed = 33
  )
  out <- do.call(get_fdelta_curves, c(list(x = x), args))
  out_renamed <- do.call(get_fdelta_curves, c(list(x = x_renamed), args))
  expect_equal(out$curves$estimate, out_renamed$curves$estimate, tolerance = 1e-12)
  expect_lt(max(abs(out$curves$reconstruction_error)), 1e-12)
})

test_that("exact and Monte Carlo tie averages agree", {
  x <- tie_example(
    reference = c(2, 2, 1), focal = c(2, 1, 2),
    mass = c(1, 10, 100)
  )
  pair <- fundecomp:::.fd_prepare_pair(
    list(series = "s", time = 1), list(series = "s", time = 2), x
  )
  exact <- fundecomp:::.fd_subset_state_variants(
    pair, "traits", x, tie_method = "average",
    max_exact_tie_perm = 10000, tie_seed = 4
  )
  monte_carlo <- fundecomp:::.fd_subset_state_variants(
    pair, "traits", x, tie_method = "average",
    n_tie_perm = 5000, max_exact_tie_perm = 1, tie_seed = 4
  )
  exact_value <- fundecomp:::.fd_eval_tie_curve(
    exact, fundecomp:::.fd_get_metric("CWM"), 0, "full", NULL, 1, 5
  )$replicates
  mc_value <- fundecomp:::.fd_eval_tie_curve(
    monte_carlo, fundecomp:::.fd_get_metric("CWM"), 0, "full", NULL, 1, 5
  )$replicates
  expect_true(exact$diagnostics$exact)
  expect_false(monte_carlo$diagnostics$exact)
  expect_equal(as.numeric(mc_value), as.numeric(exact_value), tolerance = 0.1)
  expect_gt(monte_carlo$diagnostics$n_used, 1000)
})

test_that("Shapley remains efficient and label-invariant with ties", {
  x <- tie_example(
    reference = c(10, 10, 5, 5),
    focal = c(5, 5, 10, 10),
    mass = c(1, 2, 100, 101)
  )
  x_renamed <- renamed_tie_example(x, c(2, 4, 1, 3))
  args <- list(
    reference = list(series = "s", time = 1),
    focal = list(series = "s", time = 2), metric = "biomass",
    effort_unit = "full", method = "shapley", n_perm = 120,
    tie_method = "average", seed = 7
  )
  out <- do.call(get_fdelta_curves, c(list(x = x), args))
  out_renamed <- do.call(get_fdelta_curves, c(list(x = x_renamed), args))
  expect_equal(out$curves$estimate, out_renamed$curves$estimate, tolerance = 1e-12)
  expect_equal(sum(out$curves$estimate), 990, tolerance = 1e-12)
  expect_lt(max(abs(out$curves$reconstruction_error)), 1e-12)
  expect_equal(length(out$permutations), 120L)
})

test_that("tie Monte Carlo averaging is reproducible and local to RNG state", {
  x <- tie_example(
    reference = rep(5, 6), focal = rep(5, 6),
    mass = c(1, 2, 4, 8, 16, 32)
  )
  args <- list(
    reference = list(series = "s", time = 1),
    focal = list(series = "s", time = 2), metric = "CWM",
    effort_unit = "full", tie_method = "average", n_tie_perm = 20L,
    max_exact_tie_perm = 1L, seed = 13, tie_seed = 29
  )
  set.seed(101)
  before <- .Random.seed
  first <- do.call(get_fdelta_curves, c(list(x = x), args))
  expect_identical(before, .Random.seed)
  second <- do.call(get_fdelta_curves, c(list(x = x), args))
  expect_identical(first$curves, second$curves)
  expect_equal(first$diagnostics$tie$states[[2]]$n_used, 20L)
  expect_false(first$diagnostics$tie$states[[2]]$exact)
})

test_that("Shapley reruns retain explicit public tie settings", {
  x <- tie_example(
    reference = rep(5, 4), focal = rep(5, 4), mass = c(1, 2, 4, 8)
  )
  path <- decompose_fmetric_change(
    x, list(series = "s", time = 1), list(series = "s", time = 2),
    metric = "CWM", trait = "mass", effort_unit = "full",
    n_rep = 1L, seed = 5L, tie_method = "average", n_tie_perm = 7L,
    max_exact_tie_perm = 1L, tie_seed = 41L
  )

  shapley <- shapley_decomposition(path, exact = FALSE, n_orders = 3L, seed = 8L)

  expect_identical(shapley$diagnostics$tie$n_tie_perm, 7L)
  expect_identical(shapley$diagnostics$tie$max_exact_tie_perm, 1L)
  expect_identical(shapley$diagnostics$tie$seed, 41L)
})
