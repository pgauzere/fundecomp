test_that("kNCN is rejected for longitude-latitude coordinates", {
  comm <- rbind(p1 = c(a = 2, b = 1), p2 = c(a = 1, b = 2),
                p3 = c(a = 1, b = 1))
  traits <- data.frame(species = c("a", "b"), mass = c(1, 2))
  coords <- data.frame(longitude = c(-156.5, -156.4, -156.3),
                       latitude = c(21.4, 21.5, 21.6))
  x <- make_fundecomp_input(comm = comm, traits = traits, plot_attr = coords,
                    coord_names = c("longitude", "latitude"), latlong = TRUE)

  expect_error(
    rarefy_community(x, "CWM", 1:2, effort_unit = "sample", spat_algo = "kNCN"),
    "experimental.*not supported with longitude/latitude"
  )
  expect_no_error(
    rarefy_community(x, "CWM", 1:2, effort_unit = "sample", spat_algo = "kNN")
  )
})
