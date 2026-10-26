test_that("stochastic public functions restore the caller RNG state", {
  comm <- matrix(c(3, 1, 0, 0, 2, 1), nrow = 2, byrow = TRUE,
                 dimnames = list(c("p1", "p2"), c("a", "b", "c")))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4))
  x <- make_fundecomp_input(comm = comm, traits = traits)
  st <- fundecomp:::.fd_states(x)[[1L]]
  set.seed(99)
  before <- .Random.seed
  invisible(rarefy_community(st, "CWM", 1:2, "individual", seed = 10))
  expect_identical(before, .Random.seed)
  unseeded <- rarefy_community(st, "CWM", 1:2, "individual")
  expect_identical(before, .Random.seed)
  expect_true(all(is.finite(unseeded$seed)))
  invisible(null_aggregation(comm, seed = 10, n_swap = 100))
  expect_identical(before, .Random.seed)
})

test_that("individual rarefaction rejects continuous abundances", {
  comm <- matrix(c(1.5, 0, 0.5, 2), nrow = 2, byrow = TRUE,
                 dimnames = list(c("p1", "p2"), c("a", "b")))
  traits <- data.frame(species = c("a", "b"), mass = c(1, 2))
  x <- make_fundecomp_input(comm = comm, traits = traits, count_data = FALSE,
                    time = 1:2)
  st <- fundecomp:::.fd_states(x)[[1L]]
  expect_error(rarefy_community(st, "CWM", 1, "individual"), "integer count")
  out <- get_fdelta_curves(
    x, list(series = "series1", time = 1),
    list(series = "series1", time = 2), metric = "CWM",
    effort_unit = "full", seed = 5
  )
  expect_true(all(vapply(out$states, function(z) !isTRUE(z$count_data), logical(1))))
  expect_lt(max(abs(out$curves$reconstruction_error)), 1e-12)
})

test_that("fixed-margin aggregation null preserves both margins", {
  comm <- matrix(c(4, 1, 0, 0, 2, 3, 1, 0, 2), nrow = 3, byrow = TRUE,
                 dimnames = list(paste0("p", 1:3), paste0("sp", 1:3)))
  out <- null_aggregation(comm, seed = 12, n_swap = 1000)
  expect_equal(rowSums(out), rowSums(comm))
  expect_equal(colSums(out), colSums(comm))
  expect_true(all(out >= 0 & out == floor(out)))
  expect_true(all(unlist(attr(out, "constraints")[c("sample_totals", "species_totals")])))
})

test_that("pool baseline and trait weighting are operational", {
  comm <- rbind(t1 = c(a = 3, b = 0, c = 0),
                t2 = c(a = 0, b = 1, c = 0),
                t3 = c(a = 0, b = 1, c = 3))
  traits <- data.frame(species = c("a", "b", "c"), mass = c(1, 2, 10))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = rep("s", 3), time = 1:3)
  p_ref <- make_fundecomp_pool(x, baseline = "time_reference",
                          time_reference = list(series = "s", time = 1))
  expect_equal(p_ref$sad, c(a = 1, b = 0, c = 0), tolerance = 1e-12)
  p_equal <- make_fundecomp_pool(x, trait_pool = "pooled_species")
  p_weight <- make_fundecomp_pool(x, trait_pool = "abundance_weighted")
  expect_false(isTRUE(all.equal(p_equal$trait_weights, p_weight$trait_weights)))
})

test_that("dates are ordered by value rather than input row order", {
  comm <- rbind(t3 = c(a = 1, b = 1), t1 = c(a = 1, b = 1), t2 = c(a = 1, b = 1))
  traits <- data.frame(species = c("a", "b"), mass = c(1, 2))
  x <- make_fundecomp_input(comm = comm, traits = traits, time = c(3, 1, 2))
  out <- decompose_fmetric_ts(x, metric = "CWM", effort_unit = "full")
  expect_equal(unique(out$curves$time_reference), "1")
  expect_equal(unique(out$curves$time_focal), c("2", "3"))
})

test_that("multivariate STD resampling keeps observed trait rows together", {
  comm <- rbind(t1 = c(a = 2, b = 0, c = 0),
                t2 = c(a = 0, b = 1, c = 2))
  traits <- data.frame(species = c("a", "b", "c"), t1 = c(1, 2, 10), t2 = c(10, 20, 100))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = c("s", "s"), time = c(1, 2))
  pair <- fundecomp:::.fd_prepare_pair(list(series = "s", time = 1),
                                     list(series = "s", time = 2), x)
  pair$ta_axis <- fundecomp:::.fd_resolve_ta_axis(
    "t1", pair$traits, x$trait_names, c("t1", "t2")
  )
  state <- fundecomp:::.fd_subset_state(pair, "traits", x)
  source_rows <- unique(pair$traits[, c("t1", "t2")])
  observed_rows <- unique(state$traits[state$comm[, "a"] + state$comm[, "b"] + state$comm[, "c"] > 0,
                           c("t1", "t2")])
  expect_true(all(apply(observed_rows, 1, function(z)
    any(apply(source_rows, 1, function(y) isTRUE(all.equal(as.numeric(z), as.numeric(y))))))))
})

test_that("decomposition summaries retain the effort dimension", {
  comm <- rbind(t1a = c(a = 2, b = 1), t1b = c(a = 1, b = 0),
                t2a = c(a = 1, b = 2), t2b = c(a = 0, b = 1))
  traits <- data.frame(species = c("a", "b"), mass = c(1, 2))
  x <- make_fundecomp_input(comm = comm, traits = traits, time = c(1, 1, 2, 2))
  out <- get_fdelta_curves(x, list(series = "series1", time = 1),
                           list(series = "series1", time = 2), metric = "CWM",
                           effort = 1:2, effort_unit = "sample")
  stats <- get_fdelta_stats(out)
  expect_equal(nrow(stats), 10)
  expect_true(all(stats$effort %in% 1:2))
})

test_that("intermediate counterfactual states are integer with exact margins", {
  comm <- rbind(
    r1 = c(a = 4, b = 1, c = 0), r2 = c(a = 2, b = 1, c = 0),
    f1 = c(a = 1, b = 1, c = 2), f2 = c(a = 1, b = 2, c = 4)
  )
  traits <- data.frame(species = c("a", "b", "c"), mass = c(1, 2, 5))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = rep("s", 4),
                    time = c(1, 1, 2, 2))
  pair <- fundecomp:::.fd_prepare_pair(list(series = "s", time = 1),
                                     list(series = "s", time = 2), x)
  state <- fundecomp:::.fd_subset_state(pair, c("N", "SAD"), x)
  expect_true(state$count_data)
  expect_true(all(state$comm >= 0 & state$comm == floor(state$comm)))
  expect_equal(sum(state$comm), sum(pair$focal$comm))
  expect_equal(unname(rowSums(state$comm)), unname(state$construction$row_targets))
  expect_equal(unname(colSums(state$comm)), unname(state$construction$col_targets))
})

test_that("IPF converges for sparse patterns and unequal margins", {
  pattern <- matrix(0, nrow = 12, ncol = 20)
  pattern[cbind(seq_len(12), seq_len(12))] <- 1
  dimnames(pattern) <- list(paste0("p", seq_len(12)),
                            paste0("sp", seq_len(20)))
  out <- fundecomp:::.fd_integer_counterfactual(
    pattern = pattern,
    row_weights = rep(1, 12),
    col_weights = exp(-seq(0, 7, length.out = 20)),
    N = 960L
  )
  expect_equal(sum(out), 960)
  expect_equal(as.integer(rowSums(out)),
               as.integer(attr(out, "construction")$row_targets))
  expect_equal(as.integer(colSums(out)),
               as.integer(attr(out, "construction")$col_targets))
})

test_that("individual curves use one nested draw across effort", {
  comm <- matrix(c(3, 2, 2), nrow = 1,
                 dimnames = list("p1", c("a", "b", "c")))
  traits <- data.frame(species = c("a", "b", "c"), code = c(1, 10, 100))
  x <- make_fundecomp_input(comm = comm, traits = traits)
  st <- fundecomp:::.fd_states(x)[[1L]]
  encoded_sum <- function(abundance, traits, trait = NULL, ...) {
    sum(abundance * traits$code)
  }
  out <- rarefy_community(st, encoded_sum, effort = seq_len(sum(comm)),
                          effort_unit = "individual", n_rep = 8, seed = 91)
  reps <- attr(out, "replicates")
  increments <- t(apply(reps, 1, function(z) c(z[1L], diff(z))))
  expect_true(all(increments %in% traits$code))
  expect_true(isTRUE(attr(out, "nested_sampling")))
  expect_equal(reps[, ncol(reps)], rep(sum(colSums(comm) * traits$code), 8))
})

test_that("individual decomposition works on all integerized states", {
  comm <- rbind(
    r1 = c(a = 4, b = 1, c = 0), r2 = c(a = 2, b = 1, c = 0),
    f1 = c(a = 1, b = 1, c = 2), f2 = c(a = 1, b = 2, c = 4)
  )
  traits <- data.frame(species = c("a", "b", "c"), mass = c(1, 2, 5))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = rep("s", 4),
                    time = c(1, 1, 2, 2))
  out <- get_fdelta_curves(
    x, list(series = "s", time = 1), list(series = "s", time = 2),
    metric = "CWM", effort = 1:7, effort_unit = "individual",
    n_perm = 20, seed = 44
  )
  expect_true(all(vapply(out$states, function(z)
    isTRUE(z$count_data) && all(z$comm == floor(z$comm)), logical(1))))
  constructed <- Filter(function(z) !is.null(z$construction), out$states)
  expect_true(all(vapply(constructed, function(z) {
    identical(as.integer(rowSums(z$comm)),
              as.integer(z$construction$row_targets)) &&
      identical(as.integer(colSums(z$comm)),
                as.integer(z$construction$col_targets))
  }, logical(1))))
  expect_equal(dim(out$contrast_replicates), NULL)
  expect_equal(length(out$contrast_replicates), 5)
  for (j in seq_along(out$contrast_replicates)) {
    rows <- out$curves$component == names(out$contrast_replicates)[j]
    expect_equal(out$curves$estimate[rows],
                 unname(colMeans(out$contrast_replicates[[j]])))
    expect_equal(out$curves$lower[rows],
                 unname(apply(out$contrast_replicates[[j]], 2, stats::quantile,
                              probs = 0.025, names = FALSE)))
  }
  reconstructed <- Reduce("+", out$contrast_replicates)
  expect_equal(reconstructed, out$endpoint_replicates, tolerance = 1e-12)
  expect_lt(max(abs(out$curves$reconstruction_error)), 1e-12)
})
