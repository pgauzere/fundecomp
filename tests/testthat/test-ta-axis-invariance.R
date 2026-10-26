ta_axis_example <- function(column_order = c("mass", "pace"), row_order = NULL) {
  comm <- rbind(
    r1 = c(a = 8L, b = 4L, c = 2L, d = 1L),
    r2 = c(a = 2L, b = 4L, c = 5L, d = 4L),
    f1 = c(a = 1L, b = 3L, c = 6L, d = 4L),
    f2 = c(a = 2L, b = 5L, c = 4L, d = 3L)
  )
  traits <- data.frame(
    species = c("a", "b", "c", "d"),
    mass = c(1, 2, 4, 8),
    pace = c(8, 6, 3, 1),
    climate = c(-2, -1, 1, 2),
    stringsAsFactors = FALSE
  )
  if (!is.null(row_order)) traits <- traits[row_order, , drop = FALSE]
  traits <- traits[, c("species", column_order, setdiff("climate", column_order)),
                   drop = FALSE]
  make_fundecomp_input(
    comm = comm, traits = traits, series = rep("s", 4),
    time = c(1, 1, 2, 2), sample = rownames(comm)
  )
}

ta_axis_curves <- function(x, metric = "FDis", trait = c("mass", "pace"),
                           ta_axis = "mass", method = "mobr_path",
                           n_perm = 1L) {
  get_fdelta_curves(
    x, list(series = "s", time = 1), list(series = "s", time = 2),
    metric = metric, trait = trait, ta_axis = ta_axis,
    effort_unit = "full", method = method, n_perm = n_perm, seed = 101
  )
}

test_that("univariate T-A remains backward compatible", {
  comm <- rbind(reference = c(a = 8L, b = 4L, c = 2L, d = 1L),
                focal = c(a = 1L, b = 3L, c = 6L, d = 5L))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 4, 8))
  x <- make_fundecomp_input(comm = comm, traits = traits,
                    series = rep("s", 2), time = 1:2)
  args <- list(x = x, reference = list(series = "s", time = 1),
               focal = list(series = "s", time = 2), metric = "biomass",
               effort_unit = "full", seed = 22)
  implicit <- do.call(get_fdelta_curves, args)
  explicit <- do.call(get_fdelta_curves, c(args, list(ta_axis = "mass")))
  expect_equal(implicit$curves, explicit$curves, tolerance = 1e-12)

  implicit <- do.call(get_fdelta_curves, c(args, list(method = "shapley", n_perm = 120L)))
  explicit <- do.call(get_fdelta_curves, c(args, list(ta_axis = "mass",
                                                      method = "shapley", n_perm = 120L)))
  expect_equal(implicit$curves, explicit$curves, tolerance = 1e-12)
})

test_that("explicit T-A axis makes multivariate results invariant to column order", {
  x1 <- ta_axis_example(c("mass", "pace"))
  x2 <- ta_axis_example(c("pace", "mass"))
  a <- ta_axis_curves(x1)
  b <- ta_axis_curves(x2)
  cols <- c("component", "estimate", "lower", "upper", "reconstruction_error")
  expect_equal(a$curves[cols], b$curves[cols], tolerance = 1e-12)

  a <- ta_axis_curves(x1, method = "shapley", n_perm = 120L)
  b <- ta_axis_curves(x2, method = "shapley", n_perm = 120L)
  expect_equal(a$curves[cols], b$curves[cols], tolerance = 1e-12)
})

test_that("rescaling a non-axis trait does not change T-A ordering", {
  x1 <- ta_axis_example(c("mass", "pace"))
  x2 <- ta_axis_example(c("mass", "pace"))
  x2$traits$pace <- x2$traits$pace * 100
  a <- ta_axis_curves(x1, metric = "biomass", trait = "mass")
  b <- ta_axis_curves(x2, metric = "biomass", trait = "mass")
  expect_equal(a$curves, b$curves, tolerance = 1e-12)
  expect_equal(a$ta_axis$scores, b$ta_axis$scores, tolerance = 1e-12)
})

test_that("named T-A score vectors align by species, not input row or column order", {
  x1 <- ta_axis_example(c("mass", "pace"))
  x2 <- ta_axis_example(c("pace", "mass"), row_order = c(4, 2, 1, 3))
  scores <- stats::setNames(c(a = 0.1, b = 0.4, c = 0.2, d = 0.9),
                            c("a", "b", "c", "d"))
  a <- ta_axis_curves(x1, ta_axis = scores)
  b <- ta_axis_curves(x2, ta_axis = scores)
  cols <- c("component", "estimate", "lower", "upper", "reconstruction_error")
  expect_equal(a$curves[cols], b$curves[cols], tolerance = 1e-12)
  expect_equal(unname(a$ta_axis$scores), unname(scores), tolerance = 1e-12)
})

test_that("multivariate T-A requires an explicit axis", {
  x <- ta_axis_example()
  expect_error(
    ta_axis_curves(x, ta_axis = NULL),
    "Multivariate trait vectors have no unique total ordering for T-A decomposition"
  )
  expect_error(
    ta_axis_curves(x, ta_axis = stats::setNames(c(1, 2), c("a", "b"))),
    "missing scores"
  )
  expect_error(
    ta_axis_curves(x, ta_axis = c(1, 2, 3, 4)),
    "named numeric"
  )
  expect_error(
    ta_axis_curves(x, ta_axis = c(a = 1, a = 2, b = 3, c = 4, d = 5)),
    "named numeric"
  )
  expect_error(
    ta_axis_curves(x, ta_axis = c(a = 1, b = 2, c = Inf, d = 4)),
    "finite"
  )
})

test_that("T-A ordering preserves complete multivariate rows", {
  x <- ta_axis_example(c("mass", "pace"))
  out <- ta_axis_curves(x)
  state_after_traits <- out$states[[5L]]
  got <- unique(state_after_traits$traits[colSums(state_after_traits$comm) > 0,
                                          c("mass", "pace")])
  expected <- unique(x$traits[, c("mass", "pace")])
  key <- function(z) apply(z, 1, paste, collapse = "|")
  expect_equal(sort(unname(key(got))), sort(unname(key(expected))))
})

test_that("Shapley with an explicit axis is efficient and column-order invariant", {
  x1 <- ta_axis_example(c("mass", "pace"))
  x2 <- ta_axis_example(c("pace", "mass"))
  a <- ta_axis_curves(x1, method = "shapley", n_perm = 120L)
  b <- ta_axis_curves(x2, method = "shapley", n_perm = 120L)
  expect_lt(max(abs(a$curves$reconstruction_error)), 1e-12)
  expect_equal(a$curves$estimate, b$curves$estimate, tolerance = 1e-12)
  expect_equal(sum(a$curves$estimate), sum(b$curves$estimate), tolerance = 1e-12)
  reference_value <- rarefy_community(
    a$endpoint$reference, metric = "FDis", trait = c("mass", "pace"),
    effort = 0, effort_unit = "full"
  )$value
  focal_value <- rarefy_community(
    a$endpoint$focal, metric = "FDis", trait = c("mass", "pace"),
    effort = 0, effort_unit = "full"
  )$value
  expect_equal(sum(a$curves$estimate), focal_value - reference_value,
               tolerance = 1e-12)
})

test_that("multivariate explicit-axis T-A retains abundance-tie invariance", {
  comm <- rbind(reference = c(a = 10L, b = 10L, c = 5L, d = 5L),
                focal = c(a = 5L, b = 5L, c = 10L, d = 10L))
  traits <- data.frame(
    species = c("a", "b", "c", "d"),
    mass = c(1, 2, 100, 101), pace = c(8, 7, 2, 1),
    stringsAsFactors = FALSE
  )
  x <- make_fundecomp_input(comm = comm, traits = traits,
                    series = rep("s", 2), time = 1:2)
  x_renamed <- make_fundecomp_input(
    comm = comm[, c("c", "a", "d", "b"), drop = FALSE],
    traits = traits[c(3, 1, 4, 2), c("species", "pace", "mass")],
    series = rep("s", 2), time = 1:2
  )
  args <- list(reference = list(series = "s", time = 1),
               focal = list(series = "s", time = 2), metric = "FDis",
               trait = c("mass", "pace"), ta_axis = "mass",
               effort_unit = "full", tie_method = "average", seed = 19)
  a <- do.call(get_fdelta_curves, c(list(x = x), args))
  b <- do.call(get_fdelta_curves, c(list(x = x_renamed), args))
  expect_equal(a$curves$estimate, b$curves$estimate, tolerance = 1e-12)
  expect_equal(a$curves$reconstruction_error,
               b$curves$reconstruction_error, tolerance = 1e-12)
})
