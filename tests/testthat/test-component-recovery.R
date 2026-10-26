component_estimates <- function(comm, traits, times, metric = "biomass",
                                effort = NULL, effort_unit = "full",
                                n_perm = 1L, seed = 17L,
                                method = "mobr_path") {
  x <- make_fundecomp_input(
    comm = comm, traits = traits, series = rep("recovery", nrow(comm)),
    time = times, sample = rownames(comm)
  )
  out <- get_fdelta_curves(
    x,
    reference = list(series = "recovery", time = 1),
    focal = list(series = "recovery", time = 2),
    metric = metric, effort = effort, effort_unit = effort_unit,
    method = method,
    n_perm = n_perm, seed = seed
  )
  stats::setNames(out$curves$estimate, out$curves$component)
}

expect_isolated_recovery <- function(estimates, target, tolerance = 1e-8) {
  expect_gt(abs(unname(estimates[target])), tolerance)
  expect_lt(max(abs(estimates[setdiff(names(estimates), target)])), tolerance)
}

test_that("N-only change is recovered without component leakage", {
  comm <- rbind(reference = c(8, 6, 4, 2), focal = c(16, 12, 8, 4))
  colnames(comm) <- paste0("sp", 1:4)
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4, 8))
  estimate <- component_estimates(comm, traits, times = 1:2)
  expect_isolated_recovery(estimate, "N")
})

test_that("aggregation-only change is recovered at sample effort", {
  reference <- diag(12L, 4L)
  focal <- matrix(3L, nrow = 4L, ncol = 4L)
  comm <- rbind(reference, focal)
  rownames(comm) <- paste0(rep(c("r", "f"), each = 4), seq_len(4))
  colnames(comm) <- paste0("sp", 1:4)
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4, 8))
  estimate <- component_estimates(
    comm, traits, times = rep(1:2, each = 4), metric = "richness",
    effort = 1, effort_unit = "sample", n_perm = 20
  )
  expect_isolated_recovery(estimate, "aggregation")
})

test_that("SAD-only change uses an unlabelled abundance margin", {
  comm <- rbind(reference = c(10, 6, 3, 1), focal = c(7, 6, 4, 3))
  colnames(comm) <- paste0("sp", 1:4)
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4, 8))
  estimate <- component_estimates(comm, traits, times = 1:2)
  expect_isolated_recovery(estimate, "SAD")
})

test_that("STD-only turnover preserves the SAD and pairing", {
  comm <- rbind(
    reference = c(10, 6, 3, 1, 0, 0, 0, 0),
    focal = c(0, 0, 0, 0, 10, 6, 3, 1)
  )
  colnames(comm) <- letters[1:8]
  traits <- data.frame(
    species = colnames(comm),
    mass = c(1, 2, 3, 4, 2, 4, 6, 8)
  )
  estimate <- component_estimates(comm, traits, times = 1:2)
  expect_isolated_recovery(estimate, "traits")
})

test_that("multivariate STD transport keeps complete trait rows", {
  comm <- rbind(
    reference = c(10, 6, 3, 1, 0, 0, 0, 0),
    focal = c(0, 0, 0, 0, 10, 6, 3, 1)
  )
  colnames(comm) <- letters[1:8]
  traits <- data.frame(
    species = letters[1:8],
    mass = c(1, 2, 3, 4, 2, 4, 6, 8),
    pace = c(8, 6, 7, 5, 1, 4, 2, 3)
  )
  x <- make_fundecomp_input(
    comm = comm, traits = traits, series = c("s", "s"), time = 1:2
  )
  pair <- fundecomp:::.fd_prepare_pair(
    list(series = "s", time = 1), list(series = "s", time = 2), x
  )
  pair$ta_axis <- fundecomp:::.fd_resolve_ta_axis(
    "mass", pair$traits, x$trait_names, c("mass", "pace")
  )
  state <- fundecomp:::.fd_subset_state(
    pair, c("N", "aggregation", "SAD", "traits"), x
  )
  got <- state$traits[colSums(state$comm) > 0, , drop = FALSE]
  expected <- pair$traits[colSums(pair$focal$comm) > 0, , drop = FALSE]
  row_sort <- function(z) z[do.call(order, unname(as.list(z))), , drop = FALSE]
  expect_equal(unname(as.matrix(row_sort(got))),
               unname(as.matrix(row_sort(expected))))
})

test_that("focal N, aggregation, and SAD reuse the observed focal community", {
  comm <- rbind(
    ref_1 = c(a = 6, b = 2, c = 1),
    ref_2 = c(a = 4, b = 1, c = 2),
    focal_1 = c(a = 2, b = 5, c = 3),
    focal_2 = c(a = 1, b = 3, c = 6)
  )
  traits <- data.frame(species = colnames(comm), mass = c(1, 3, 8))
  x <- make_fundecomp_input(
    comm = comm, traits = traits, series = rep("s", 4), time = c(1, 1, 2, 2)
  )
  pair <- fundecomp:::.fd_prepare_pair(
    list(series = "s", time = 1), list(series = "s", time = 2), x
  )

  state <- fundecomp:::.fd_subset_state(
    pair, c("N", "aggregation", "SAD"), x
  )

  expect_identical(state$construction$method, "exact_focal_community")
  expect_identical(state$comm, pair$focal$comm)
})

test_that("Trait-Abundance-only reversal is assigned to TA", {
  abundance <- c(30, 20, 15, 10, 8, 7, 6, 4)
  comm <- rbind(reference = abundance, focal = rev(abundance))
  colnames(comm) <- letters[1:8]
  traits <- data.frame(
    species = colnames(comm), mass = 10^seq(0, 2, length.out = 8)
  )
  estimate <- component_estimates(comm, traits, times = 1:2)
  expect_isolated_recovery(estimate, "TA")

  shapley <- component_estimates(
    comm, traits, times = 1:2, method = "shapley", n_perm = 120
  )
  expect_isolated_recovery(shapley, "TA")
})
