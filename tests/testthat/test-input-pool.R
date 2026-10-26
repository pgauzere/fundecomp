test_that("input validation records and drops incomplete traits", {
  comm <- matrix(c(2, 1, 0, 1, 0, 2), nrow = 2, byrow = TRUE,
                 dimnames = list(c("p1", "p2"), c("sp1", "sp2", "sp3")))
  traits <- data.frame(species = c("sp1", "sp2"), mass = c(1, NA))
  x <- make_fundecomp_input(comm = comm, traits = traits, time = c(1, 2))
  expect_equal(x$diagnostics$dropped_species, c("sp2", "sp3"))
  expect_equal(colnames(x$comm), "sp1")
  expect_error(make_fundecomp_input(comm = comm, traits = traits, trait_na = "error"),
               "without complete traits")
})

test_that("dataset pool uses validated defaults", {
  comm <- matrix(c(2, 1, 0, 1, 0, 2), nrow = 2, byrow = TRUE,
                 dimnames = list(c("p1", "p2"), c("sp1", "sp2", "sp3")))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 3))
  x <- make_fundecomp_input(comm = comm, traits = traits, time = c(1, 2))
  p <- make_fundecomp_pool(x)
  expect_equal(p$sad_pool, "pooled_empirical")
  expect_equal(p$trait_pool, "pooled_species")
  expect_equal(p$species, colnames(comm))
  expect_equal(sum(p$sad), 1)
})

test_that("leave-one-time-out pool excludes only the selected SAD contribution", {
  comm <- matrix(c(2, 1, 0, 1, 0, 2), nrow = 2, byrow = TRUE,
                 dimnames = list(c("p1", "p2"), c("sp1", "sp2", "sp3")))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 3))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = c("s", "s"), time = c(1, 2))
  p <- make_fundecomp_pool(x, exclude_time = 2)
  expect_equal(p$baseline, "dataset_pool_leave_one_time_out")
  expect_equal(p$diagnostics$excluded_time, 2)
  expect_equal(p$species, colnames(comm))
  expect_equal(unname(p$sad), c(2 / 3, 1 / 3, 0), tolerance = 1e-12)
  expect_identical(p$diagnostics$sad_role, "descriptive_pool_summary_only")
})

test_that("descriptive pools cannot be passed as decomposition baselines", {
  comm <- rbind(
    r1 = c(a = 4, b = 1, c = 0), r2 = c(a = 2, b = 1, c = 0),
    f1 = c(a = 1, b = 1, c = 2), f2 = c(a = 1, b = 2, c = 4)
  )
  traits <- data.frame(species = c("a", "b", "c"), mass = c(1, 2, 5))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = rep("s", 4),
                    time = c(1, 1, 2, 2))
  p_empirical <- make_fundecomp_pool(x, sad_pool = "pooled_empirical")
  p_uniform <- make_fundecomp_pool(x, sad_pool = "uniform")
  args <- list(x = x, reference = list(series = "s", time = 1),
               focal = list(series = "s", time = 2), metric = "biomass",
               effort_unit = "full", method = "mobr_path", seed = 21)
  expect_false(isTRUE(all.equal(p_empirical$sad, p_uniform$sad)))
  expect_error(
    do.call(get_fdelta_curves, c(args, list(pool = p_empirical))),
    "Regional-pool SAD/STD baselines are not implemented"
  )
  expect_error(
    do.call(decompose_fmetric_change, c(args, list(pool = p_uniform))),
    "Regional-pool SAD/STD baselines are not implemented"
  )
})

test_that("input rejects invalid counts", {
  comm <- matrix(c(1, -1), nrow = 1, dimnames = list("p1", c("sp1", "sp2")))
  traits <- data.frame(species = c("sp1", "sp2"), mass = c(1, 2))
  expect_error(make_fundecomp_input(comm = comm, traits = traits), "non-negative")
  expect_error(make_fundecomp_input(comm = matrix(c(1.5, 0), nrow = 1,
                                           dimnames = list("p1", c("sp1", "sp2"))),
                             traits = traits), "integer")
})

test_that("long input permits the same site at multiple dates", {
  long <- data.frame(
    series = rep("s", 4), time = c(1, 1, 2, 2),
    sample = rep("site1", 4), species = rep(c("sp1", "sp2"), 2),
    abundance = c(2, 1, 1, 3)
  )
  traits <- data.frame(species = c("sp1", "sp2"), mass = c(1, 2))
  x <- make_fundecomp_input(long = long, traits = traits)
  expect_equal(nrow(x$comm), 2)
  expect_equal(length(unique(x$meta$time)), 2)
})
