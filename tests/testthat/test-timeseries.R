test_that("time series defaults to the first date", {
  comm <- matrix(c(2, 1, 1, 2, 0, 3, 3, 0, 1), nrow = 3, byrow = TRUE,
                 dimnames = list(paste0("p", 1:3), paste0("sp", 1:3)))
  traits <- data.frame(species = colnames(comm), mass = c(1, 2, 5))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = rep("s", 3), time = c(2001, 2003, 2008))
  out <- decompose_fmetric_ts(x, metric = "biomass", effort_unit = "full")
  expect_s3_class(out, "fundecomp_timeseries")
  expect_true(all(out$curves$time_reference == "2001"))
  expect_true(all(out$curves$time_focal %in% c("2003", "2008")))
  expect_true(all(out$curves$delta_t %in% c(2, 7)))
})
