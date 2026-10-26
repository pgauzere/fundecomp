unequal_richness_input <- function(s_ref = 8L, s_focal = 12L,
                                   tied = FALSE, multivariate = FALSE,
                                   spatial = FALSE) {
  species <- paste0("sp", seq_len(max(s_ref, s_focal)))
  ref <- focal <- stats::setNames(rep.int(0L, length(species)), species)
  ref[seq_len(s_ref)] <- if (tied) rep(c(8L, 4L), length.out = s_ref) else seq(s_ref, 1L)
  focal[seq_len(s_focal)] <- if (tied) rep(c(8L, 4L), length.out = s_focal) else seq(s_focal, 1L)
  if (spatial) {
    comm <- rbind(ref_1 = ref * c(rep(1L, ceiling(length(ref) / 2)), rep(0L, floor(length(ref) / 2)))[seq_along(ref)],
                  ref_2 = ref - ref * c(rep(1L, ceiling(length(ref) / 2)), rep(0L, floor(length(ref) / 2)))[seq_along(ref)],
                  focal_1 = focal * c(rep(1L, ceiling(length(focal) / 2)), rep(0L, floor(length(focal) / 2)))[seq_along(focal)],
                  focal_2 = focal - focal * c(rep(1L, ceiling(length(focal) / 2)), rep(0L, floor(length(focal) / 2)))[seq_along(focal)])
    time <- c(1, 1, 2, 2)
  } else {
    comm <- rbind(reference = ref, focal = focal)
    time <- 1:2
  }
  traits <- data.frame(species = species, mass = seq_along(species),
                       stringsAsFactors = FALSE)
  if (multivariate) {
    traits$pace <- sin(seq_along(species))
    traits$habitat <- seq_along(species)^2
  }
  make_fundecomp_input(comm = comm, traits = traits, series = rep("s", nrow(comm)),
               time = time, sample = rownames(comm))
}

unequal_curves <- function(x, metric = "biomass", method = "mobr_path",
                           effort_unit = "full", effort = NULL, ...) {
  get_fdelta_curves(
    x, list(series = "s", time = 1), list(series = "s", time = 2),
    metric = metric,
    trait = if (metric %in% c("FDis", "FRic")) c("mass", "pace") else
      if (metric %in% c("biomass", "CWM", "variance")) "mass" else NULL,
    ta_axis = if ("pace" %in% names(x$traits)) "mass" else NULL,
    method = method, effort = effort, effort_unit = effort_unit,
    n_perm = if (method == "shapley") 120L else 1L, seed = 41, ...
  )
}

test_that("unequal richness reconstructs path and Shapley endpoints", {
  for (sizes in list(c(8L, 12L), c(12L, 8L))) {
    x <- unequal_richness_input(sizes[1], sizes[2])
    for (metric in c("biomass", "CWM", "richness")) {
      path <- unequal_curves(x, metric)
      shapley <- unequal_curves(x, metric, method = "shapley")
      expect_lt(max(abs(path$curves$reconstruction_error)), 1e-12)
      expect_lt(max(abs(shapley$curves$reconstruction_error)), 1e-12)
      expect_equal(sum(shapley$curves$estimate),
                   sum(shapley$curves$estimate) - shapley$curves$reconstruction_error[1],
                   tolerance = 1e-12)
    }
  }
  x <- unequal_richness_input(8L, 12L, multivariate = TRUE)
  expect_lt(max(abs(unequal_curves(x, "FDis")$curves$reconstruction_error)), 1e-12)
})

test_that("unequal-richness states are invariant to arbitrary labels with abundance ties", {
  x <- unequal_richness_input(8L, 12L, tied = TRUE)
  perm <- c(12, 3, 9, 1, 7, 4, 11, 2, 10, 5, 8, 6)
  renamed_comm <- x$comm[, perm, drop = FALSE]
  colnames(renamed_comm) <- paste0("renamed_", seq_len(ncol(renamed_comm)))
  renamed_traits <- data.frame(species = colnames(renamed_comm),
                               mass = x$traits$mass[perm])
  renamed <- make_fundecomp_input(comm = renamed_comm, traits = renamed_traits,
                           series = x$meta$series, time = x$meta$time,
                           sample = x$meta$sample)
  args <- list(metric = "biomass", tie_method = "average")
  a <- do.call(unequal_curves, c(list(x = x), args))
  b <- do.call(unequal_curves, c(list(x = renamed), args))
  expect_equal(a$curves$estimate, b$curves$estimate, tolerance = 1e-12)
  expect_lt(max(abs(a$curves$reconstruction_error)), 1e-12)
})

test_that("reverse unequal-richness comparison reverses endpoints and Shapley efficiency", {
  x <- unequal_richness_input(8L, 12L)
  forward <- unequal_curves(x, "biomass", method = "shapley")
  backward <- get_fdelta_curves(
    x, list(series = "s", time = 2), list(series = "s", time = 1),
    metric = "biomass", effort_unit = "full", method = "shapley",
    trait = "mass", n_perm = 120L, seed = 41
  )
  expect_equal(sum(forward$curves$estimate), -sum(backward$curves$estimate),
               tolerance = 1e-12)
  expect_lt(max(abs(backward$curves$reconstruction_error)), 1e-12)
})

test_that("multivariate unequal-richness transport preserves full rows and T-A axis semantics", {
  x <- unequal_richness_input(8L, 12L, multivariate = TRUE)
  a <- unequal_curves(x, "FDis", method = "shapley")
  swapped <- x
  swapped$traits <- swapped$traits[, c("habitat", "pace", "mass"), drop = FALSE]
  b <- unequal_curves(swapped, "FDis", method = "shapley")
  rescaled <- x
  rescaled$traits$pace <- rescaled$traits$pace * 1000
  # Rescaling a non-axis trait can change FDis itself; biomass based on the
  # axis trait isolates the requested invariant T-A ordering diagnostic.
  c <- unequal_curves(rescaled, "biomass", method = "shapley")
  d <- unequal_curves(x, "biomass", method = "shapley")
  expect_equal(a$curves$estimate, b$curves$estimate, tolerance = 1e-12)
  expect_equal(d$curves$estimate, c$curves$estimate, tolerance = 1e-12)
  expect_lt(max(abs(a$curves$reconstruction_error)), 1e-12)
})

test_that("unequal richness remains defined at partial spatial effort", {
  x <- unequal_richness_input(8L, 12L, tied = TRUE, spatial = TRUE)
  out <- unequal_curves(x, "CWM", effort_unit = "sample", effort = 1:2,
                        tie_method = "average", n_tie_perm = 20L)
  expect_true(all(is.finite(out$curves$estimate)))
  expect_lt(max(abs(out$curves$reconstruction_error)), 1e-12)
})
