plot_example <- function(itv = FALSE) {
  comm <- rbind(a = c(sp1 = 3, sp2 = 1, sp3 = 0),
                b = c(sp1 = 1, sp2 = 2, sp3 = 0),
                c = c(sp1 = 1, sp2 = 1, sp3 = 2),
                d = c(sp1 = 0, sp2 = 2, sp3 = 2))
  traits <- data.frame(species = colnames(comm), size = c(1, 2, 5))
  x <- make_fundecomp_input(comm = comm, traits = traits,
                    series = rep("s", 4), time = c(1, 1, 2, 2))
  get_fdelta_curves(x, list(series = "s", time = 1), list(series = "s", time = 2),
                    metric = "trait_sum", effort_unit = "full", seed = 4)
}

test_that("plot types return ggplot objects and preserve decomposition values", {
  out <- plot_example()
  original <- out$curves
  expect_s3_class(plot(out), "ggplot")
  expect_s3_class(plot(out, type = "curves"), "ggplot")
  expect_s3_class(plot(out, type = "contributions"), "ggplot")
  expect_s3_class(plot(out, type = "contributions", effort = 10), "ggplot")
  expect_s3_class(plot(out, type = "reconstruction"), "ggplot")
  expect_no_warning(ggplot2::ggplot_build(plot(out, type = "reconstruction")))
  expect_s3_class(plot(out) + ggplot2::geom_point(), "ggplot")
  expect_identical(out$curves, original)
  expect_error(plot(out, type = "unknown"), "arg")
})

test_that("full-effort contribution plot is not labelled as effort zero", {
  out <- plot_example()
  expect_identical(plot(out, type = "contributions")$labels$title,
                   "Full sampling effort")
})

test_that("contribution plot includes ITV when present", {
  out <- plot_example()
  out$curves <- rbind(out$curves, transform(out$curves[1, ], component = "ITV",
                                            component_label = "ITV", estimate = 0,
                                            lower = 0, upper = 0))
  out$path <- c(out$path, "ITV")
  expect_s3_class(plot(out, type = "contributions"), "ggplot")
})
