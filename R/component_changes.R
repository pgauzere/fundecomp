# ggplot2 uses tidy evaluation for column names in plot mappings.
utils::globalVariables(c("state", "N", "relative_abundance", "abundance",
                         "value", ".data", "site_rank", "relative_site_abundance",
                         "species", "density", "rank_for_plot"))

#' Extract observed ecological component states
#'
#' Extracts endpoint descriptions from a `fundecomp_input` or retained
#' `fundecomp_result`. These are observed endpoint states, not hybrid states
#' created during decomposition. Rank-based SAD and TA tables expose the
#' package's rank-slot representation; aggregation is returned as the observed
#' sample-by-species abundance matrix in long form.
#'
#' @param x A `fundecomp_input` or `fundecomp_result`.
#' @param reference,focal Endpoint identifiers, each a list with `series` and
#'   `time` (or `state`). Required for an input object and default to the
#'   retained endpoints for a result object.
#' @param components Components to extract: `N`, `SAD`, `aggregation`, `STD`,
#'   `TA`, and `ITV` when individual trait observations are available.
#' @param ta_axis Optional trait column name or named numeric species scores
#'   used to order trait rows for TA. Multivariate traits require this argument.
#' @return A named list containing endpoint metadata and one data frame per
#'   selected component. Component rows identify their endpoint in `state`.
#' @details STD contains one unweighted species trait row per species present
#'   at the endpoint. TA pairs trait rows with abundance ranks ordered by
#'   `ta_axis`; this describes the observed association and does not reproduce
#'   tie-averaged hybrid mappings. ITV contains individual trait residuals
#'   around endpoint-specific species centroids. Missing species slots are not
#'   invented. Aggregation is shown as observed sample allocation; the
#'   reconstruction transition additionally rank-transports profiles and
#'   fits margins, so it is not itself a conventional aggregation index.
#' @examples
#' # See the package vignette for a complete two-endpoint example.
#' @export
extract_component_changes <- function(x, reference = NULL, focal = NULL,
                                      components = c("N", "SAD", "aggregation", "STD", "TA"),
                                      ta_axis = NULL) {
  if (inherits(x, "fundecomp_result")) {
    if (is.null(x$.input)) .fd_abort("The result does not retain its input object.")
    if (is.null(reference)) reference <- x$reference
    if (is.null(focal)) focal <- x$focal
    x <- x$.input
  }
  if (!inherits(x, "fundecomp_input")) .fd_abort("`x` must be a fundecomp_input or fundecomp_result.")
  if (is.null(reference) || is.null(focal)) .fd_abort("Supply `reference` and `focal` endpoints.")
  if (is.null(reference$time)) reference$time <- reference$state
  if (is.null(focal$time)) focal$time <- focal$state
  available <- c("N", "SAD", "aggregation", "STD", "TA", if (isTRUE(x$has_itv)) "ITV")
  components[components == "traits"] <- "STD"
  components[components == "Aggregation"] <- "aggregation"
  if (!length(components) || anyDuplicated(components) || any(!components %in% available)) {
    .fd_abort("`components` must be selected from: ", paste(available, collapse = ", "), ".")
  }
  # Build endpoint states with the same input state constructor used by the
  # decomposition. The union support is useful for aligned species tables.
  states <- list(reference = .fd_state(x, reference$series, reference$time),
                 focal = .fd_state(x, focal$series, focal$time))
  species_union <- union(states$reference$species, states$focal$species)
  traits_by_state <- lapply(states, `[[`, "traits")
  axis_traits <- x$traits
  axis <- .fd_resolve_ta_axis(ta_axis, axis_traits, x$trait_names,
                              if (is.null(ta_axis) && length(x$trait_names) == 1L) x$trait_names else NULL)
  answer <- list(metadata = list(reference = reference, focal = focal,
                                 data_type = "observed endpoint states",
                                 components = components, ta_axis = axis$name))
  if ("N" %in% components) {
    answer$N <- data.frame(state = names(states), N = vapply(states, `[[`, numeric(1), "N"),
                           stringsAsFactors = FALSE)
    answer$N$change <- answer$N$N[2] - answer$N$N[1]
    answer$N$relative_change_percent <- if (answer$N$N[1] != 0) 100 * answer$N$change / answer$N$N[1] else NA_real_
  }
  if ("SAD" %in% components) {
    rows <- list()
    for (state_name in names(states)) {
      abundance <- colSums(states[[state_name]]$comm)
      abundance <- abundance[abundance > 0]
      abundance <- sort(as.numeric(abundance), decreasing = TRUE)
      rows[[state_name]] <- data.frame(state = state_name, rank = seq_along(abundance),
                                       abundance = abundance,
                                       relative_abundance = abundance / sum(abundance))
    }
    answer$SAD <- do.call(rbind, rows)
  }
  if ("aggregation" %in% components) {
    rows <- list()
    for (state_name in names(states)) {
      # Matrix row names in temporal inputs may encode series and time. The
      # metadata sample ID is the actual sampling-unit identity required for
      # comparing site allocation between observed endpoints.
      state_community <- states[[state_name]]$comm
      rownames(state_community) <- as.character(states[[state_name]]$meta$sample)
      z <- as.data.frame(as.table(state_community), stringsAsFactors = FALSE)
      names(z) <- c("sample", "species", "abundance")
      z$state <- state_name
      rows[[state_name]] <- z[c("state", "sample", "species", "abundance")]
    }
    answer$aggregation <- do.call(rbind, rows)
  }
  if ("STD" %in% components || "TA" %in% components) {
    trait_rows <- list(); ta_rows <- list()
    for (state_name in names(states)) {
      present <- colnames(states[[state_name]]$comm)[colSums(states[[state_name]]$comm) > 0]
      td <- states[[state_name]]$traits[present, , drop = FALSE]
      if ("STD" %in% components) {
        z <- data.frame(state = state_name, species = present, td, check.names = FALSE)
        trait_rows[[state_name]] <- z
      }
      if ("TA" %in% components) {
        abund <- colSums(states[[state_name]]$comm)[present]
        rank_order <- .fd_rank_order(abund)
        axis_order <- .fd_ta_axis_order(rank_order, axis$scores)
        z <- data.frame(state = state_name, rank = seq_along(rank_order),
                        species = rank_order, abundance = as.numeric(abund[rank_order]),
                        axis_rank = match(rank_order, axis_order),
                        td[rank_order, , drop = FALSE], check.names = FALSE)
        ta_rows[[state_name]] <- z
      }
    }
    if ("STD" %in% components) answer$STD <- do.call(rbind, trait_rows)
    if ("TA" %in% components) answer$TA <- do.call(rbind, ta_rows)
  }
  if ("ITV" %in% components) {
    rows <- list()
    for (state_name in names(states)) for (sp in names(states[[state_name]]$itv)) {
      residual <- as.data.frame(states[[state_name]]$itv[[sp]])
      residual$state <- state_name; residual$species <- sp
      residual$observation <- seq_len(nrow(residual))
      rows[[length(rows) + 1L]] <- residual[c("state", "species", "observation", x$trait_names)]
    }
    answer$ITV <- do.call(rbind, rows)
  }
  class(answer) <- c("fundecomp_component_changes", "list")
  answer
}

#' Plot observed ecological component changes
#'
#' Creates one ggplot per selected component from extracted component data or
#' directly from a compatible input/result object. With `combine = TRUE`,
#' the returned list also contains a patchwork multipanel plot.
#' @param x Component data, a `fundecomp_input`, or a `fundecomp_result`.
#' @param components Optional component selection when `x` is an input/result.
#' @param reference,focal Endpoint identifiers when `x` is an input object.
#' @param traits Optional trait columns to display for STD and TA.
#' @param colours Two colours named `reference` and `focal`.
#' @param state_labels Display labels named `reference` and `focal`. Defaults
#'   are `Initial state` and `Final state`; customise these for non-temporal comparisons.
#' @param sad_style SAD display: grouped bars or rank-abundance lines. When
#'   omitted, lines are used for more than 30 ranks and bars for smaller SADs.
#' @param aggregation_style Spatial rank profiles (default), ranked-site
#'   heatmap, observed sample-allocation bars, or observed site-identity heatmap.
#' @param log_rank Use a logarithmic x-axis for rank-based SAD curves. This is
#'   most useful with `sad_style = "line"`.
#' @param ta_style Trait values by abundance rank as points (default) or bars.
#' @param log_abundance Use a logarithmic relative-abundance axis for SAD.
#'   When omitted, it is enabled with the automatic line display for long SADs.
#' @param ta_axis Optional scalar axis for multivariate TA extraction.
#' @param combine Also return a patchwork multipanel plot as `combined`.
#' @return A named list of individual ggplot objects and, optionally, a combined plot.
#' @details SAD rank curves use numeric ranks; `log_rank` and `log_abundance`
#'   change only the display axes. Zero-abundance species are absent from an
#'   observed SAD, so the logarithmic relative-abundance axis has no zero values.
#'   The spatial heatmap shows each present species' observed abundance share
#'   at each sampling unit. Its site-identity ordering is shared across states;
#'   absent species have undefined shares. Rank-profile plots and ranked-site
#'   heatmaps remove site identity and show concentration only. The ranked
#'   heatmap retains each species profile without averaging across species;
#'   endpoint colours identify state, while tile opacity encodes its share on
#'   a common 0--1 scale.
#'   Neither spatial endpoint plot is the rank transport and margin-fitting
#'   operation used to construct hybrids.
#'   STD density curves use equal species weights and a bandwidth estimated
#'   from pooled endpoint values for each trait. TA points retain individual
#'   trait values; tied abundance ranks share their midrank x position rather
#'   than receiving a visual species order within the tie.
#' @export
plot_component_changes <- function(x, components = NULL, reference = NULL, focal = NULL,
                                   traits = NULL,
                                   colours = c(reference = "#0072B2", focal = "#D55E00"),
                                   state_labels = c(reference = "Initial state", focal = "Final state"),
                                   sad_style = NULL,
                                   aggregation_style = c("rank", "sample", "heatmap", "rank_heatmap"),
                                   log_rank = FALSE, ta_style = c("point", "bar"),
                                   log_abundance = NULL, ta_axis = NULL,
                                   combine = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) .fd_abort("Install ggplot2 to use plot_component_changes().")
  if (!inherits(x, "fundecomp_component_changes")) {
    if (is.null(components)) components <- c("N", "SAD", "aggregation", "STD", "TA", if (inherits(x, "fundecomp_input") && x$has_itv) "ITV")
    x <- extract_component_changes(x, reference, focal, components, ta_axis)
  }
  if (length(colours) != 2L) .fd_abort("`colours` must contain reference and focal colours.")
  if (is.null(names(colours))) names(colours) <- c("reference", "focal")
  if (!all(c("reference", "focal") %in% names(colours))) {
    .fd_abort("`colours` must be named `reference` and `focal`.")
  }
  if (is.null(names(state_labels)) || !all(c("reference", "focal") %in% names(state_labels))) {
    .fd_abort("`state_labels` must be named `reference` and `focal`.")
  }
  # Long SADs become hard to read as grouped bars. Use lines and a log
  # relative-abundance axis automatically beyond 30 ranks. This changes only
  # the visual display, not the extracted component data.
  maximum_sad_ranks <- if (!is.null(x$SAD) && nrow(x$SAD)) max(x$SAD$rank) else 0L
  if (is.null(sad_style)) {
    if (maximum_sad_ranks > 30L) {
      sad_style <- "line"
    } else {
      sad_style <- "bar"
    }
  }
  if (is.null(log_abundance)) {
    log_abundance <- maximum_sad_ranks > 30L && sad_style == "line"
  }
  sad_style <- match.arg(sad_style, c("bar", "line"))
  aggregation_style <- match.arg(aggregation_style)
  ta_style <- match.arg(ta_style)
  if (isTRUE(log_rank) && sad_style != "line") {
    .fd_abort("Use `sad_style = \"line\"` when requesting `log_rank = TRUE`.")
  }
  state_levels <- unname(state_labels[c("reference", "focal")])
  state_colours <- stats::setNames(unname(colours[c("reference", "focal")]), state_levels)
  plots <- list()
  common_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))

  if (!is.null(x$N)) {
    n_data <- x$N
    n_data$state <- factor(state_labels[n_data$state], levels = state_levels)
    initial_n <- n_data$N[1L]
    final_n <- n_data$N[2L]
    relative_change <- if (initial_n == 0) NA_real_ else 100 * (final_n - initial_n) / initial_n
    change_text <- if (is.na(relative_change)) "Relative change undefined (initial N = 0)" else
      paste0("Change: ", if (relative_change > 0) "+" else "", format(round(relative_change, 1), nsmall = 1), "%")
    plots$N <- ggplot2::ggplot(n_data, ggplot2::aes(state, N, fill = state)) +
      ggplot2::geom_col(width = 0.72, show.legend = FALSE) + ggplot2::scale_fill_manual(values = state_colours) +
      ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.08))) +
      ggplot2::labs(x = NULL, y = "Total abundance (N)", title = "Total abundance", subtitle = change_text) + common_theme
  }
  if (!is.null(x$SAD)) {
    sad_data <- x$SAD
    sad_data$state <- factor(state_labels[sad_data$state], levels = state_levels)
    # Keep ranks numeric so optional logarithmic visual scales work correctly.
    plots$SAD <- ggplot2::ggplot(sad_data, ggplot2::aes(rank, relative_abundance))
    if (sad_style == "bar") {
      # A missing rank stays missing when richness differs between endpoints.
      plots$SAD <- plots$SAD + ggplot2::geom_col(ggplot2::aes(fill = state),
        position = ggplot2::position_dodge(width = 0.82), width = 0.76)
    } else {
      # One continuous curve per state keeps long rank-abundance series legible.
      plots$SAD <- plots$SAD + ggplot2::geom_line(ggplot2::aes(colour = state, group = state), linewidth = 0.85)
    }
    if (sad_style == "bar") {
      plots$SAD <- plots$SAD + ggplot2::scale_fill_manual(values = state_colours)
    } else {
      plots$SAD <- plots$SAD + ggplot2::scale_colour_manual(values = state_colours, name = NULL)
    }
    plots$SAD <- plots$SAD +
      ggplot2::labs(x = "Abundance rank", y = "Relative abundance", title = "Species abundance distribution") + common_theme
    if (isTRUE(log_rank)) {
      plots$SAD <- plots$SAD + ggplot2::scale_x_log10(
        breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500))
    } else if (sad_style == "bar") {
      sad_rank_breaks <- pretty(range(sad_data$rank), n = 6)
      sad_rank_breaks <- sort(unique(c(1, sad_rank_breaks[sad_rank_breaks >= 1])))
      # Leave room for the first bar's half-width while keeping the first
      # displayed tick at rank one rather than suggesting a rank-zero species.
      plots$SAD <- plots$SAD + ggplot2::scale_x_continuous(
        limits = c(0.5, NA), breaks = sad_rank_breaks,
        expand = ggplot2::expansion(mult = c(0, 0.02)))
    } else {
      sad_rank_breaks <- pretty(range(sad_data$rank), n = 6)
      sad_rank_breaks <- sort(unique(c(1, sad_rank_breaks[sad_rank_breaks >= 1])))
      plots$SAD <- plots$SAD + ggplot2::scale_x_continuous(
        limits = c(1, NA), breaks = sad_rank_breaks,
        expand = ggplot2::expansion(mult = c(0, 0.02)))
    }
    if (isTRUE(log_abundance)) plots$SAD <- plots$SAD + ggplot2::scale_y_log10()
  }
  if (!is.null(x$aggregation)) {
    if (aggregation_style == "sample") {
      sample_data <- x$aggregation
      sample_data$state <- factor(state_labels[sample_data$state], levels = state_levels)
      plots$aggregation <- ggplot2::ggplot(sample_data, ggplot2::aes(sample, abundance, fill = state)) +
        ggplot2::geom_col(position = "dodge") + ggplot2::facet_wrap(~ species, scales = "free_y") +
        ggplot2::scale_fill_manual(values = state_colours) +
        ggplot2::labs(x = "Sampling unit", y = "Observed abundance", title = "Spatial allocation") + common_theme
    } else if (aggregation_style == "heatmap") {
      # A site-identity heatmap preserves which monitored site held each
      # species' share of individuals. Species order is pooled abundance rank;
      # site order is common alphanumeric sample ID order in both years.
      heatmap_data <- x$aggregation
      heatmap_data$state <- factor(state_labels[heatmap_data$state], levels = state_levels)
      profile_key <- paste(heatmap_data$state, heatmap_data$species, sep = "@@")
      heatmap_data$relative_site_abundance <- NA_real_
      for (key in unique(profile_key)) {
        selected <- which(profile_key == key)
        species_total <- sum(heatmap_data$abundance[selected])
        if (species_total > 0) {
          heatmap_data$relative_site_abundance[selected] <-
            heatmap_data$abundance[selected] / species_total
        }
      }
      pooled_species_abundance <- stats::aggregate(
        abundance ~ species, data = heatmap_data, FUN = sum)
      species_order <- pooled_species_abundance$species[
        order(-pooled_species_abundance$abundance, pooled_species_abundance$species)]
      site_order <- sort(unique(as.character(heatmap_data$sample)))
      heatmap_data$species <- factor(heatmap_data$species, levels = rev(species_order))
      heatmap_data$sample <- factor(heatmap_data$sample, levels = site_order)
      site_label_breaks <- site_order[unique(round(seq(1, length(site_order), length.out = min(8L, length(site_order)))))]
      plots$aggregation <- ggplot2::ggplot(heatmap_data,
        ggplot2::aes(sample, species, fill = relative_site_abundance)) +
        ggplot2::geom_tile(na.rm = TRUE) +
        ggplot2::facet_grid(. ~ state) +
        ggplot2::scale_fill_viridis_c(trans = "sqrt", limits = c(0, 1), breaks = c(0, 0.5, 1),
          na.value = "grey92", name = "Species share (square-root scale)",
          guide = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5,
            barwidth = grid::unit(5, "cm"), barheight = grid::unit(0.3, "cm"))) +
        ggplot2::scale_x_discrete(breaks = site_label_breaks) +
        ggplot2::labs(x = "Monitoring site (ordered by site ID)", y = "Species (ordered by pooled abundance)",
                      title = "Spatial allocation",
                      subtitle = "Each cell is the fraction of a species' individuals at a site") +
        common_theme + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5.5),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          panel.background = ggplot2::element_rect(fill = "grey92", colour = NA),
          strip.text = ggplot2::element_text(face = "bold"))
    } else if (aggregation_style == "rank_heatmap") {
      # Calculate each species' share at each sampling unit, then rank sites
      # independently within that species and endpoint. Zero-abundance sites
      # remain in each profile. Site identities are no longer shown, so tied
      # values can be ordered arbitrarily without changing the displayed data.
      ranked_heatmap_data <- x$aggregation
      ranked_heatmap_data$state <- factor(state_labels[ranked_heatmap_data$state],
                                           levels = state_levels)
      ranked_heatmap_data$relative_site_abundance <- NA_real_
      ranked_heatmap_data$site_rank <- NA_integer_
      profile_key <- paste(ranked_heatmap_data$state, ranked_heatmap_data$species, sep = "@@")
      for (key in unique(profile_key)) {
        selected <- which(profile_key == key)
        species_total <- sum(ranked_heatmap_data$abundance[selected])
        if (species_total > 0) {
          proportions <- ranked_heatmap_data$abundance[selected] / species_total
          ranked_order <- order(-proportions, seq_along(proportions))
          ranked_heatmap_data$relative_site_abundance[selected] <- proportions
          ranked_heatmap_data$site_rank[selected[ranked_order]] <- seq_along(selected)
        }
      }
      pooled_species_abundance <- stats::aggregate(
        abundance ~ species, data = ranked_heatmap_data, FUN = sum)
      species_order <- pooled_species_abundance$species[
        order(-pooled_species_abundance$abundance, pooled_species_abundance$species)]
      ranked_heatmap_data$species <- factor(ranked_heatmap_data$species,
                                             levels = rev(species_order))
      maximum_site_rank <- if (any(!is.na(ranked_heatmap_data$site_rank))) {
        max(ranked_heatmap_data$site_rank, na.rm = TRUE)
      } else {
        1L
      }
      site_rank_breaks <- pretty(c(1, maximum_site_rank), n = 5)
      site_rank_breaks <- sort(unique(c(1, site_rank_breaks[site_rank_breaks >= 1])))
      plots$aggregation <- ggplot2::ggplot(ranked_heatmap_data,
        ggplot2::aes(site_rank, species, fill = state, alpha = relative_site_abundance)) +
        ggplot2::geom_tile(na.rm = TRUE) +
        ggplot2::facet_grid(. ~ state) +
        ggplot2::scale_fill_manual(values = state_colours, name = NULL) +
        ggplot2::scale_alpha_continuous(trans = "sqrt", limits = c(0, 1),
          range = c(0.08, 1), breaks = c(0, 0.25, 1), name = "Species share",
          guide = ggplot2::guide_legend(title.position = "top", title.hjust = 0.5,
            nrow = 1, byrow = TRUE)) +
        ggplot2::scale_x_continuous(limits = c(0.5, NA), breaks = site_rank_breaks,
          expand = ggplot2::expansion(mult = c(0, 0.01))) +
        ggplot2::labs(x = "Site rank within species (largest allocation first)",
          y = "Species (ordered by pooled abundance)", title = "Spatial allocation",
          subtitle = "Site identity is discarded; species share is shown by tile opacity") +
        common_theme + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5.5),
          panel.grid = ggplot2::element_blank(),
          panel.background = ggplot2::element_rect(fill = "grey92", colour = NA),
          strip.text = ggplot2::element_text(face = "bold"))
    } else {
      # Rank sites separately within each species and endpoint. Zero-abundance
      # sites remain in profiles for species present at that endpoint. A
      # species absent from an endpoint has no defined relative allocation.
      profile_rows <- list()
      for (endpoint in c("reference", "focal")) {
        endpoint_data <- x$aggregation[x$aggregation$state == endpoint, , drop = FALSE]
        for (species_name in unique(endpoint_data$species)) {
          species_data <- endpoint_data[endpoint_data$species == species_name, , drop = FALSE]
          species_total <- sum(species_data$abundance)
          if (species_total <= 0) next
          species_data$relative_site_abundance <- species_data$abundance / species_total
          site_order <- order(-species_data$relative_site_abundance, seq_len(nrow(species_data)))
          species_data <- species_data[site_order, , drop = FALSE]
          species_data$site_rank <- seq_len(nrow(species_data))
          species_data$state <- state_labels[[endpoint]]
          profile_rows[[length(profile_rows) + 1L]] <-
            species_data[c("state", "species", "site_rank", "relative_site_abundance")]
        }
      }
      if (length(profile_rows)) {
        profile_data <- do.call(rbind, profile_rows)
        profile_data$state <- factor(profile_data$state, levels = state_levels)
        spatial_plot <- ggplot2::ggplot(profile_data,
          ggplot2::aes(site_rank, relative_site_abundance, colour = state,
                       group = interaction(state, species)))
        if (length(unique(profile_data$species)) <= 12L) {
          spatial_plot <- spatial_plot + ggplot2::geom_line(linewidth = 0.7) + ggplot2::facet_wrap(~ species)
        } else {
          # For many species, overlay each transparent profile. No across-
          # species mean or other summary is introduced.
          spatial_plot <- spatial_plot + ggplot2::geom_line(linewidth = 0.5, alpha = 0.12)
        }
        plots$aggregation <- spatial_plot + ggplot2::scale_colour_manual(values = state_colours) +
          ggplot2::labs(x = "Site abundance rank", y = "Proportion of species individuals",
                        title = "Spatial allocation",
                        subtitle = "Observed within-species concentration across sampling units") + common_theme
      }
    }
  }
  selected_traits <- traits
  if (is.null(selected_traits) && !is.null(x$STD)) selected_traits <- setdiff(names(x$STD), c("state", "species"))
  if (is.null(selected_traits) && !is.null(x$TA)) selected_traits <- setdiff(names(x$TA), c("state", "rank", "species", "abundance", "axis_rank"))
  if (!is.null(x$STD) && length(selected_traits)) {
    # Base reshape keeps the optional plotting dependency limited to ggplot2.
    long <- stats::reshape(x$STD[c("state", "species", selected_traits)],
                           varying = selected_traits, v.names = "value",
                           timevar = "trait", times = selected_traits,
                           direction = "long")
    long$state <- factor(state_labels[long$state], levels = state_levels)
    # Calculate each trait's densities on a shared grid with one bandwidth
    # estimated from both endpoint samples pooled together. This keeps the
    # two endpoint curves directly comparable while leaving species weights
    # equal. Ranks or abundance do not enter this calculation.
    density_rows <- list()
    for (trait_name in selected_traits) {
      trait_data <- long[long$trait == trait_name, , drop = FALSE]
      pooled_values <- trait_data$value[is.finite(trait_data$value)]
      if (length(pooled_values) < 2L) next
      common_bandwidth <- stats::bw.nrd0(pooled_values)
      if (!is.finite(common_bandwidth) || common_bandwidth <= 0) next
      density_grid <- seq(min(pooled_values) - 3 * common_bandwidth,
                          max(pooled_values) + 3 * common_bandwidth,
                          length.out = 256L)
      for (state_name in state_levels) {
        state_values <- trait_data$value[trait_data$state == state_name]
        state_values <- state_values[is.finite(state_values)]
        if (!length(state_values)) next
        # An explicit Gaussian kernel also handles a single species value.
        density_value <- vapply(density_grid, function(grid_value) {
          sum(stats::dnorm((grid_value - state_values) / common_bandwidth)) /
            (length(state_values) * common_bandwidth)
        }, numeric(1))
        density_rows[[length(density_rows) + 1L]] <- data.frame(
          state = state_name, trait = trait_name, value = density_grid,
          density = density_value
        )
      }
    }
    density_data <- if (length(density_rows)) do.call(rbind, density_rows) else
      data.frame(state = character(), trait = character(), value = numeric(), density = numeric())
    plots$STD <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = density_data,
        ggplot2::aes(value, ymin = 0, ymax = density, fill = state, group = state),
        alpha = 0.12, show.legend = FALSE) +
      ggplot2::geom_line(data = density_data,
        ggplot2::aes(value, density, colour = state, group = state), linewidth = 0.8) +
      ggplot2::geom_rug(data = long, ggplot2::aes(value, colour = state),
        sides = "b", alpha = 0.65, show.legend = FALSE) +
      ggplot2::facet_wrap(~ trait, scales = "free") +
      ggplot2::scale_colour_manual(values = state_colours, name = NULL) + ggplot2::scale_fill_manual(values = state_colours) +
      ggplot2::guides(fill = "none") +
      ggplot2::labs(x = "Species trait value (unweighted)", y = "Density", title = "Species trait distribution") + common_theme
  }
  if (!is.null(x$TA) && length(selected_traits)) {
    trait <- selected_traits[1]
    ta_data <- x$TA
    ta_data$state <- factor(state_labels[ta_data$state], levels = state_levels)
    # Tied abundances occupy the same midrank on this descriptive plot. Every
    # species trait value remains a separate point; only the x position is
    # shared, avoiding a false species ordering within abundance ties.
    ta_data$rank_for_plot <- stats::ave(ta_data$abundance, ta_data$state,
      FUN = function(state_abundance) rank(-state_abundance, ties.method = "average"))
    ta_rank_breaks <- pretty(c(1, max(ta_data$rank_for_plot)), n = 6)
    ta_rank_breaks <- sort(unique(c(1, ta_rank_breaks[ta_rank_breaks >= 1])))
    plots$TA <- ggplot2::ggplot(ta_data, ggplot2::aes(rank_for_plot, .data[[trait]]))
    if (ta_style == "bar") {
      plots$TA <- plots$TA + ggplot2::geom_col(ggplot2::aes(fill = state), width = 0.72)
    } else {
      plots$TA <- plots$TA + ggplot2::geom_point(ggplot2::aes(colour = state), size = 1.25, alpha = 0.45)
    }
    if (ta_style == "bar") {
      plots$TA <- plots$TA + ggplot2::scale_fill_manual(values = state_colours)
    } else {
      plots$TA <- plots$TA + ggplot2::scale_colour_manual(values = state_colours, name = NULL)
    }
    plots$TA <- plots$TA +
      ggplot2::labs(x = "Abundance rank (1 = most abundant; ties share midrank)", y = paste("Trait value:", trait),
                    title = "Trait-abundance association",
                    subtitle = "Observed species trait values; tied abundance ranks share x positions") + common_theme
    plots$TA <- plots$TA + ggplot2::scale_x_continuous(
      limits = c(if (ta_style == "bar") 0.5 else 1, NA), breaks = ta_rank_breaks,
      expand = ggplot2::expansion(mult = c(0.01, 0.02)))
  }
  if (!is.null(x$ITV) && length(selected_traits)) {
    trait <- selected_traits[1]
    itv_data <- x$ITV
    itv_data$state <- factor(state_labels[itv_data$state], levels = state_levels)
    plots$ITV <- ggplot2::ggplot(itv_data, ggplot2::aes(.data[[trait]], colour = state, fill = state)) +
      ggplot2::geom_density(alpha = 0.18, na.rm = TRUE) + ggplot2::facet_wrap(~ species) +
      ggplot2::scale_colour_manual(values = state_colours) + ggplot2::scale_fill_manual(values = state_colours) +
      ggplot2::guides(fill = "none") +
      ggplot2::labs(x = "Within-species trait deviation", y = "Density", title = "Intraspecific trait variation") + common_theme
  }
  if (isTRUE(combine) && length(plots)) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      .fd_abort("Install the optional `patchwork` package to create a combined figure.")
    }
    # Keep one discrete endpoint legend in the combined figure. Prefer TA,
    # because its points make the state colours easiest to identify. The
    # numeric spatial heatmap colour bar is a separate guide and is retained.
    combined_plots <- plots
    state_legend_component <- if (!is.null(combined_plots$TA)) {
      "TA"
    } else if (!is.null(combined_plots$SAD)) {
      "SAD"
    } else if (!is.null(combined_plots$STD)) {
      "STD"
    } else if (!is.null(combined_plots$ITV)) {
      "ITV"
    } else if (!is.null(combined_plots$aggregation)) {
      "aggregation"
    } else {
      NULL
    }
    if (!is.null(combined_plots$SAD) && state_legend_component != "SAD") {
      combined_plots$SAD <- combined_plots$SAD + ggplot2::theme(legend.position = "none")
    }
    if (!is.null(combined_plots$STD) && state_legend_component != "STD") {
      combined_plots$STD <- combined_plots$STD + ggplot2::theme(legend.position = "none")
    }
    if (!is.null(combined_plots$ITV) && state_legend_component != "ITV") {
      combined_plots$ITV <- combined_plots$ITV + ggplot2::theme(legend.position = "none")
    }
    if (!is.null(combined_plots$aggregation) && state_legend_component != "aggregation") {
      if (aggregation_style == "rank" || aggregation_style == "sample") {
        combined_plots$aggregation <- combined_plots$aggregation + ggplot2::theme(legend.position = "none")
      } else if (aggregation_style == "rank_heatmap") {
        # Keep the continuous species-share key but use the single endpoint
        # legend from TA rather than a second blue/orange heatmap key.
        combined_plots$aggregation <- combined_plots$aggregation + ggplot2::guides(fill = "none")
      }
    }
    if (!is.null(combined_plots$aggregation) && aggregation_style == "rank_heatmap" &&
        length(unique(combined_plots$aggregation$data$species)) > 30L) {
      # Species names remain in the individual plot. In a compact combined
      # figure, hiding dozens of overlapping names makes the heatmap easier to
      # read while its species rows and shared ordering remain unchanged.
      combined_plots$aggregation <- combined_plots$aggregation +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())
    }
    # patchwork assembles all selected plots and collects their remaining
    # endpoint key together with any separate continuous spatial scale.
    plots$combined <- patchwork::wrap_plots(combined_plots, ncol = 2, guides = "collect")
  }
  plots
}
