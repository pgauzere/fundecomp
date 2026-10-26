# Optional extraction and visualisation of observed component changes
# Run this script line by line in RStudio from any working directory.

# ------------------------------------------------------------
# Load the package and prepare a small reproducible example
# ------------------------------------------------------------
library(fundecomp)

# Rows are sampling units and columns are species. The first two rows are
# reference samples and the last two are focal samples. These labels are
# chronological here, so the default display labels Initial and Final apply.
community <- rbind(
  reference_1 = c(sp1 = 5, sp2 = 2, sp3 = 1),
  reference_2 = c(sp1 = 2, sp2 = 3, sp3 = 1),
  focal_1 = c(sp1 = 1, sp2 = 3, sp3 = 5),
  focal_2 = c(sp1 = 1, sp2 = 2, sp3 = 4)
)

species_traits <- data.frame(
  species = c("sp1", "sp2", "sp3"),
  mass = c(1, 3, 8)
)

component_input <- make_fundecomp_input(
  comm = community,
  traits = species_traits,
  series = rep("example_site", 4),
  state = c("initial", "initial", "final", "final")
)

# ------------------------------------------------------------
# Extract endpoint component states without a metric calculation
# ------------------------------------------------------------
component_data <- extract_component_changes(
  component_input,
  reference = list(series = "example_site", state = "initial"),
  focal = list(series = "example_site", state = "final"),
  components = c("N", "SAD", "aggregation", "STD", "TA")
)

# Inspect simple tables and then inspect all data with str().
component_data$N
head(component_data$SAD)
head(component_data$aggregation)
component_data$STD
component_data$TA
str(component_data)

# ------------------------------------------------------------
# Plot individual components and a combined multipanel figure
# ------------------------------------------------------------
component_plots <- plot_component_changes(
  component_data,
  traits = "mass",
  colours = c(reference = "#2C7FB8", focal = "#D95F0E"),
  combine = TRUE
)

# N compares observed totals and reports percentage change.
component_plots$N
# SAD bars compare relative abundance by species rank, independently of N.
component_plots$SAD
# Spatial allocation ranks sites independently within each species and state.
component_plots$aggregation
# STD compares unweighted species trait values; rug marks show each species.
component_plots$STD
# TA shows the trait value paired with each observed abundance rank.
component_plots$TA
# The combined patchwork figure prints directly in RStudio.
component_plots$combined

# ------------------------------------------------------------
# Save figures and export a component table
# ------------------------------------------------------------
output_directory <- file.path(getwd(), "component_changes_results")
dir.create(output_directory, showWarnings = FALSE)

ggplot2::ggsave(
  filename = file.path(output_directory, "sad_rank_abundance.png"),
  plot = component_plots$SAD,
  width = 6, height = 4, dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_directory, "trait_distribution.png"),
  plot = component_plots$STD,
  width = 6, height = 4, dpi = 300
)

# Save the combined patchwork figure directly.
ggplot2::ggsave(
  filename = file.path(output_directory, "component_changes_combined.pdf"),
  plot = component_plots$combined,
  width = 10, height = 8
)

utils::write.csv(
  component_data$SAD,
  file = file.path(output_directory, "sad_component_data.csv"),
  row.names = FALSE
)

# The returned plots and data remain ordinary R objects for further analysis.
