test_that("spatial richness kNN agrees with mobr sSBR", {
  comm <- matrix(c(4, 0, 1, 2, 1, 0, 0, 3, 2, 1, 0, 1, 0, 2, 3, 0),
                nrow = 4, byrow = TRUE,
                dimnames = list(paste0("p", 1:4), paste0("sp", 1:4)))
  traits <- data.frame(species = colnames(comm), trait = 1:4)
  coords <- data.frame(x = c(0, 2, 0, 2), y = c(0, 0, 3, 3))
  x <- make_fundecomp_input(comm = comm, traits = traits, plot_attr = coords,
                    coord_names = c("x", "y"))
  st <- fundecomp:::.fd_states(x)[[1L]]
  f <- rarefy_community(st, "richness", 1:4, "sample", spat_algo = "kNN")$value
  m <- mobr::rarefaction(comm, method = "sSBR", effort = 1:4,
                         coords = coords, latlong = FALSE, spat_algo = "kNN")
  expect_equal(f, as.numeric(m), tolerance = 1e-12)
})

test_that("longitude-latitude kNN uses the mobr spherical distance", {
  comm <- matrix(c(3, 0, 1, 1, 0, 2, 0, 2, 1), nrow = 3, byrow = TRUE,
                 dimnames = list(paste0("p", 1:3), paste0("sp", 1:3)))
  traits <- data.frame(species = colnames(comm), trait = 1:3)
  coords <- data.frame(longitude = c(5.0, 5.2, 5.05),
                       latitude = c(46.0, 46.0, 46.2))
  x <- make_fundecomp_input(comm = comm, traits = traits, plot_attr = coords,
                    coord_names = c("longitude", "latitude"), latlong = TRUE)
  st <- fundecomp:::.fd_states(x)[[1L]]
  f <- fundecomp::rarefy_community(st, "richness", 1:3, "sample", spat_algo = "kNN")$value
  m <- mobr::rarefaction(comm, method = "sSBR", effort = 1:3,
                         coords = coords, latlong = TRUE, spat_algo = "kNN")
  expect_equal(f, as.numeric(m), tolerance = 1e-12)
})

test_that("public decomposition labels distinguish STD and Trait-Abundance", {
  comm <- matrix(c(2, 1, 0, 1, 0, 2), nrow = 2, byrow = TRUE,
                 dimnames = list(c("p1", "p2"), c("sp1", "sp2", "sp3")))
  traits <- data.frame(species = colnames(comm), trait = c(1, 2, 4))
  x <- make_fundecomp_input(comm = comm, traits = traits, series = rep("s", 2), time = c(1, 2))
  out <- get_fdelta_curves(x, list(series = "s", time = 1),
                           list(series = "s", time = 2), metric = "FDis")
  expect_true(all(c("STD", "Trait-Abundance") %in% unique(out$curves$component_label)))
})

test_that("Trait-Abundance null preserves SAD and multivariate STD", {
  comm <- matrix(c(8, 2, 1, 0, 1, 4), nrow = 2, byrow = TRUE,
                 dimnames = list(c("p1", "p2"), c("sp1", "sp2", "sp3")))
  traits <- data.frame(t1 = c(1, 2, 5), t2 = c(10, 20, 50),
                       row.names = colnames(comm))
  x <- make_fundecomp_input(comm = comm,
                    traits = data.frame(species = rownames(traits), traits),
                    species_col = "species")
  st <- fundecomp:::.fd_states(x)[[1L]]
  nul <- break_trait_abundance_link(st, permutation = c(2, 1, 3))
  expect_true(all(unlist(nul$constraints)))
  expect_equal(colSums(st$comm), colSums(nul$state$comm))
  expect_equal(sort(st$traits$t1), sort(nul$state$traits$t1))
  expect_equal(sort(st$traits$t2), sort(nul$state$traits$t2))
  expect_false(identical(st$traits, nul$state$traits))
  before <- rarefy_community(st, "CWM", 0, "full", trait = "t1")$value
  after <- rarefy_community(nul$state, "CWM", 0, "full", trait = "t1")$value
  expect_false(isTRUE(all.equal(before, after)))
})
