test_that("component extraction describes observed endpoint states", {
  comm <- rbind(r1 = c(a = 4, b = 1, c = 0), r2 = c(a = 2, b = 1, c = 0),
                f1 = c(a = 1, b = 2, c = 3), f2 = c(a = 1, b = 1, c = 2))
  traits <- data.frame(species = c("a", "b", "c"), mass = c(1, 3, 8))
  x <- make_fundecomp_input(comm, traits = traits,
                            series = rep("s", 4), state = c("r", "r", "f", "f"))
  d <- extract_component_changes(x, list(series = "s", state = "r"),
                                 list(series = "s", state = "f"),
                                 components = c("N", "SAD", "aggregation", "STD", "TA"))
  expect_equal(d$N$N, c(8, 10))
  expect_equal(d$SAD$relative_abundance[d$SAD$state == "reference"], c(0.75, 0.25))
  expect_equal(d$STD$state, c(rep("reference", 2), rep("focal", 3)))
  expect_true(all(c("rank", "axis_rank", "mass") %in% names(d$TA)))
  expect_equal(sum(d$aggregation$abundance[d$aggregation$state == "focal"]), 10)
  expect_error(extract_component_changes(x, components = "ITV"), "Supply `reference` and `focal`")
})

test_that("component extraction retains optional ITV residual data", {
  comm <- rbind(r = c(a = 2, b = 1), f = c(a = 1, b = 2))
  traits <- data.frame(species = c("a", "b"), mass = c(1, 4))
  observations <- data.frame(series = rep("s", 4), time = c(1, 1, 2, 2),
                             species = c("a", "a", "a", "b"), mass = c(0, 2, 2, 3))
  x <- make_fundecomp_input(comm, traits = traits, series = c("s", "s"), time = c(1, 2),
                            trait_observations = observations, itv_na = "species_mean")
  d <- extract_component_changes(x, list(series = "s", time = 1),
                                 list(series = "s", time = 2), components = "ITV")
  expect_true("ITV" %in% names(d))
  expect_true(all(c("state", "species", "observation", "mass") %in% names(d$ITV)))
  ta_data <- extract_component_changes(x, list(series = "s", time = 1),
                                       list(series = "s", time = 2), components = "TA")
  expect_true(all(c("mass", "rank", "axis_rank") %in% names(ta_data$TA)))
})

test_that("plotting selected components returns ggplot objects", {
  skip_if_not_installed("ggplot2")
  comm <- rbind(r = c(a = 2, b = 1), f = c(a = 1, b = 2))
  traits <- data.frame(species = c("a", "b"), mass = c(1, 4))
  x <- make_fundecomp_input(comm, traits = traits, series = c("s", "s"), state = c(1, 2))
  d <- extract_component_changes(x, list(series = "s", state = 1),
                                 list(series = "s", state = 2), components = c("N", "SAD"))
  plots <- plot_component_changes(d)
  expect_named(plots, c("N", "SAD"))
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_s3_class(plots$N$layers[[1]]$geom, "GeomCol")
  expect_s3_class(plots$SAD$layers[[1]]$geom, "GeomCol")
  expect_equal(levels(plots$N$data$state), c("Initial state", "Final state"))
  n_bars <- ggplot2::ggplot_build(plots$N)$data[[1]]
  expect_equal(n_bars$fill, c("#0072B2", "#D55E00"))
  expect_equal(plots$N$scales$get_scales("y")$limits[1], 0)
})

test_that("spatial profiles rank sites within each species and endpoint", {
  skip_if_not_installed("ggplot2")
  comm <- rbind(r1 = c(a = 9, b = 0), r2 = c(a = 1, b = 2), r3 = c(a = 0, b = 8),
                f1 = c(a = 0, b = 9), f2 = c(a = 1, b = 1), f3 = c(a = 9, b = 0))
  traits <- data.frame(species = c("a", "b"), mass = c(1, 2))
  x <- make_fundecomp_input(comm, traits = traits, series = rep("s", 6),
                            state = c(rep("r", 3), rep("f", 3)))
  d <- extract_component_changes(x, list(series = "s", state = "r"),
                                 list(series = "s", state = "f"), components = "aggregation")
  plots <- plot_component_changes(d)
  profiles <- plots$aggregation$data
  initial_a <- profiles[profiles$state == "Initial state" & profiles$species == "a", ]
  initial_b <- profiles[profiles$state == "Initial state" & profiles$species == "b", ]
  expect_equal(initial_a$relative_site_abundance, c(0.9, 0.1, 0))
  expect_equal(initial_b$relative_site_abundance, c(0.8, 0.2, 0))
  expect_equal(nrow(initial_a), 3L) # the zero-abundance sampling unit is retained
  expect_s3_class(plots$aggregation$layers[[1]]$geom, "GeomLine")
  expect_match(plots$aggregation$labels$x, "Site abundance rank")
  sample_plot <- plot_component_changes(d, aggregation_style = "sample")$aggregation
  expect_s3_class(sample_plot$layers[[1]]$geom, "GeomCol")

  heatmap_plot <- plot_component_changes(d, aggregation_style = "heatmap")$aggregation
  expect_s3_class(heatmap_plot$layers[[1]]$geom, "GeomTile")
  expect_equal(levels(heatmap_plot$data$state), c("Initial state", "Final state"))
  expect_setequal(as.character(heatmap_plot$data$sample), unique(d$aggregation$sample))

  ranked_heatmap <- plot_component_changes(d, aggregation_style = "rank_heatmap")$aggregation
  expect_s3_class(ranked_heatmap$layers[[1]]$geom, "GeomTile")
  ranked_initial_a <- ranked_heatmap$data[
    ranked_heatmap$data$state == "Initial state" & ranked_heatmap$data$species == "a", ]
  expect_equal(ranked_initial_a$site_rank, c(1L, 2L, 3L))
  expect_equal(ranked_initial_a$relative_site_abundance, c(0.9, 0.1, 0))
  expect_equal(nrow(ranked_initial_a), 3L)
  expect_equal(ranked_heatmap$scales$get_scales("alpha")$limits, c(0, 1))
  ranked_heatmap_built <- ggplot2::ggplot_build(ranked_heatmap)$data[[1]]
  expect_true(any(ranked_heatmap_built$alpha == 0.08)) # zero-share tiles remain visible
  expect_equal(unname(ranked_heatmap$scales$get_scales("fill")$palette(2)),
               c("#0072B2", "#D55E00"))
})

test_that("SAD log axes use positive numeric ranks and preserve endpoint curves", {
  skip_if_not_installed("ggplot2")
  comm <- rbind(r1 = c(a = 4, b = 1, c = 0), r2 = c(a = 2, b = 1, c = 0),
                f1 = c(a = 1, b = 2, c = 3), f2 = c(a = 1, b = 1, c = 2))
  traits <- data.frame(species = c("a", "b", "c"), mass = c(1, 3, 8))
  x <- make_fundecomp_input(comm, traits = traits, series = rep("s", 4),
                            state = c("r", "r", "f", "f"))
  d <- extract_component_changes(x, list(series = "s", state = "r"),
                                 list(series = "s", state = "f"), components = "SAD")
  plot <- plot_component_changes(d, sad_style = "line", log_rank = TRUE,
                                 log_abundance = TRUE)$SAD
  expect_true(is.numeric(plot$data$rank))
  expect_true(all(plot$data$rank >= 1))
  expect_true(all(plot$data$relative_abundance > 0))
  expect_equal(length(unique(ggplot2::ggplot_build(plot)$data[[1]]$group)), 2L)
  expect_error(plot_component_changes(d, log_rank = TRUE), "sad_style = \\\"line\\\"")
})

test_that("long SADs default to linear-rank lines with a log abundance axis", {
  skip_if_not_installed("ggplot2")
  species_names <- paste0("sp", seq_len(35))
  comm <- matrix(rep(seq_len(35), times = 4), nrow = 4,
                 dimnames = list(c("r1", "r2", "f1", "f2"), species_names))
  traits <- data.frame(species = species_names, trait = seq_len(35))
  x <- make_fundecomp_input(comm, traits = traits, series = rep("s", 4),
                            state = c("r", "r", "f", "f"))
  d <- extract_component_changes(x, list(series = "s", state = "r"),
                                 list(series = "s", state = "f"), components = "SAD")
  plot <- plot_component_changes(d)$SAD
  expect_s3_class(plot$layers[[1]]$geom, "GeomLine")
  expect_equal(plot$scales$get_scales("x")$trans$name, "identity")
  expect_equal(plot$scales$get_scales("y")$trans$name, "log-10")
  expect_equal(length(unique(ggplot2::ggplot_build(plot)$data[[1]]$group)), 2L)
})

test_that("STD overlays unweighted endpoint trait values and TA bars retain values", {
  skip_if_not_installed("ggplot2")
  comm <- rbind(r1 = c(a = 4, b = 1), r2 = c(a = 1, b = 4),
                f1 = c(a = 1, b = 4), f2 = c(a = 4, b = 1))
  traits <- data.frame(species = c("a", "b"), trait = c(-2, 5))
  x <- make_fundecomp_input(comm, traits = traits, series = rep("s", 4),
                            state = c("r", "r", "f", "f"))
  d <- extract_component_changes(x, list(series = "s", state = "r"),
                                 list(series = "s", state = "f"), components = c("STD", "TA"))
  plots <- plot_component_changes(d, traits = "trait")
  rug_data <- plots$STD$layers[[3]]$data
  expect_equal(as.integer(table(rug_data$state)), c(2L, 2L))
  expect_equal(sort(rug_data$value), c(-2, -2, 5, 5))
  expect_s3_class(plots$STD$layers[[1]]$geom, "GeomRibbon")
  expect_s3_class(plots$STD$layers[[2]]$geom, "GeomLine")
  expect_s3_class(plots$STD$layers[[3]]$geom, "GeomRug")
  expect_s3_class(plots$TA$layers[[1]]$geom, "GeomPoint")
  expect_equal(sort(plots$TA$data$trait), c(-2, -2, 5, 5))
  expect_match(plots$TA$labels$y, "trait")
  expect_equal(plots$TA$data$rank_for_plot, c(1.5, 1.5, 1.5, 1.5))
  bar_plot <- plot_component_changes(d, traits = "trait", ta_style = "bar")$TA
  expect_s3_class(bar_plot$layers[[1]]$geom, "GeomCol")
})

test_that("combined patchwork includes every requested panel and can be saved", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  comm <- rbind(r1 = c(a = 3, b = 1), r2 = c(a = 1, b = 3),
                f1 = c(a = 1, b = 3), f2 = c(a = 3, b = 1))
  traits <- data.frame(species = c("a", "b"), mass = c(1, 4))
  x <- make_fundecomp_input(comm, traits = traits, series = rep("s", 4),
                            state = c("r", "r", "f", "f"))
  d <- extract_component_changes(x, list(series = "s", state = "r"),
                                 list(series = "s", state = "f"),
                                 components = c("N", "SAD", "aggregation", "STD", "TA"))
  plots <- plot_component_changes(d, traits = "mass", combine = TRUE)
  individual_names <- c("N", "SAD", "aggregation", "STD", "TA")
  expect_true(all(individual_names %in% names(plots)))
  expect_s3_class(plots$combined, "patchwork")
  # patchwork stores the first panel as the primary plot and the remaining
  # panels in patches$plots. Counting both verifies all five are assembled.
  expect_equal(length(plots$combined$patches$plots) + 1L, length(individual_names))
  combined_titles <- c(plots$combined$labels$title,
    vapply(plots$combined$patches$plots, function(panel) panel$labels$title, character(1)))
  expect_setequal(combined_titles, c("Total abundance", "Species abundance distribution",
                                     "Spatial allocation", "Species trait distribution",
                                     "Trait-abundance association"))
  expect_no_error(ggplot2::ggplotGrob(plots$combined))
  output_file <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(output_file, plot = plots$combined, width = 8, height = 7)
  expect_gt(file.info(output_file)$size, 0)
})
