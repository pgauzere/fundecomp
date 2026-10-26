# ggplot2 mappings use these data columns.
utils::globalVariables(c("component", "effort", "estimate", "lower", "upper",
                         "component_label", "series", "value", "y", "label"))

#' Plot an `fundecomp_result` decomposition
#'
#' The default `"curves"` view shows one panel per component across sampling
#' effort. `"contributions"` shows signed contributions at one effort (the
#' largest available by default), and `"reconstruction"` compares the
#' conditional contrast and observed change across effort. Returns a ggplot
#' object that can be further customised with ggplot2 layers.
#'
#' @param x An `fundecomp_result` object.
#' @param type Plot type: component curves, contributions at one effort, or
#'   reconstruction diagnostic. Defaults to `"curves"` for compatibility.
#' @param effort Optional effort for `type = "contributions"`. The closest
#'   available effort is used.
#' @param ... Additional arguments retained for compatibility; currently
#'   unused.
#' @examples
#' # Given a decomposition result `out`:
#' # plot(out)
#' # plot(out, type = "contributions")
#' # plot(out, type = "reconstruction")
#' @importFrom ggplot2 ggplot
#' @export
plot.fundecomp_result <- function(x, type = c("curves", "contributions", "reconstruction"),
                          effort = NULL, ...) {
  if (!inherits(x, "fundecomp_result")) .fd_abort("`x` must be an `fundecomp_result` object.")
  type <- match.arg(type)
  z <- x$curves
  required <- c("component", "effort", "estimate")
  if (!is.data.frame(z) || !all(required %in% names(z)) || !nrow(z)) {
    .fd_abort("`x$curves` must contain component, effort, and estimate values.")
  }
  comp_order <- unique(as.character(z$component))
  if (!is.null(x$path)) comp_order <- c(intersect(x$path, comp_order), setdiff(comp_order, x$path))
  labels <- if ("component_label" %in% names(z)) {
    stats::setNames(as.character(z$component_label[match(comp_order, z$component)]), comp_order)
  } else stats::setNames(comp_order, comp_order)
  z$component <- factor(z$component, levels = comp_order)
  z$component_label <- factor(labels[as.character(z$component)], levels = unname(labels))
  common_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))
  component_palette <- c(
    N = "#0072B2", aggregation = "#E69F00", SAD = "#009E73",
    traits = "#CC79A7", STD = "#CC79A7", ITV = "#56B4E9", TA = "#D55E00"
  )
  component_colours <- unname(component_palette[comp_order])
  names(component_colours) <- comp_order
  missing_colours <- is.na(component_colours)
  if (any(missing_colours)) {
    component_colours[missing_colours] <- grDevices::hcl.colors(sum(missing_colours), "Dark 3")
  }

  if (type == "curves") {
    p <- ggplot2::ggplot(z, ggplot2::aes(effort, estimate, colour = component,
                                         group = component))
    if (all(c("lower", "upper") %in% names(z))) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper,
                                                  fill = component), alpha = 0.16,
                                     colour = NA, show.legend = FALSE, na.rm = TRUE)
    }
    p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
      ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
      ggplot2::facet_wrap(~ component_label, ncol = 1, scales = "free_y") +
      ggplot2::scale_colour_manual(values = component_colours, guide = "none") +
      ggplot2::scale_fill_manual(values = component_colours) +
      ggplot2::labs(x = "Sampling effort", y = "Contribution") + common_theme
    return(p)
  }

  if (type == "contributions") {
    efforts <- sort(unique(z$effort))
    if (is.null(effort)) selected <- max(efforts) else {
      if (length(effort) != 1L || !is.finite(effort)) .fd_abort("`effort` must be one finite value.")
      selected <- efforts[which.min(abs(efforts - effort))]
    }
    d <- z[z$effort == selected, , drop = FALSE]
    d <- d[match(comp_order, as.character(d$component)), , drop = FALSE]
    if (!nrow(d) || anyNA(d$estimate)) .fd_abort("Contributions are unavailable at the selected effort.")
    d$component_label <- factor(labels[as.character(d$component)], levels = unname(labels))
    p <- ggplot2::ggplot(d, ggplot2::aes(estimate, component_label, colour = component)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4)
    if (all(c("lower", "upper") %in% names(d))) {
      p <- p + ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper),
                                       orientation = "y", width = 0.18, linewidth = 0.6,
                                       na.rm = TRUE)
    }
    effort_title <- if (identical(x$effort_unit, "full")) {
      "Full sampling effort"
    } else {
      paste(tools::toTitleCase(x$effort_unit), "effort:", selected)
    }
    return(p + ggplot2::geom_point(size = 2.2) +
      ggplot2::scale_colour_manual(values = component_colours, guide = "none") +
      ggplot2::labs(x = "Signed contribution", y = NULL,
                    title = effort_title) + common_theme)
  }

  if (!"reconstruction_error" %in% names(z)) {
    .fd_abort("Reconstruction plotting requires `reconstruction_error` in `x$curves`.")
  }
  by_effort <- split(z, z$effort)
  efforts <- as.numeric(names(by_effort))
  errors <- vapply(by_effort, function(d) d$reconstruction_error[1L], numeric(1))
  reconstructed <- vapply(by_effort, function(d) sum(d$estimate), numeric(1))
  contrast_label <- if (isTRUE(x$full_decomposition)) "Observed" else "Partial counterfactual"
  diagnostic <- data.frame(effort = efforts, contrast = reconstructed - errors,
                           `Reconstructed change` = reconstructed,
                           check.names = FALSE)
  names(diagnostic)[2L] <- contrast_label
  if (!isTRUE(x$full_decomposition)) {
    diagnostic$`Observed total change` <- vapply(by_effort,
      function(d) d$observed_total_change[1L], numeric(1))
  }
  diagnostic <- utils::stack(diagnostic[-1L])
  names(diagnostic) <- c("value", "series")
  diagnostic$effort <- rep(efforts, times = if (isTRUE(x$full_decomposition)) 2L else 3L)
  line_levels <- if (isTRUE(x$full_decomposition)) c("Observed", "Reconstructed change") else
    c("Partial counterfactual", "Reconstructed change", "Observed total change")
  diagnostic$series <- factor(diagnostic$series, levels = line_levels)
  line_colours <- stats::setNames(c("#0072B2", "#D55E00", "#4D4D4D")[seq_along(line_levels)], line_levels)
  p <- ggplot2::ggplot(diagnostic, ggplot2::aes(effort, value, colour = series, linetype = series)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4)
  if (length(unique(diagnostic$effort)) > 1L) {
    p <- p + ggplot2::geom_line(linewidth = 0.8)
  }
  p + ggplot2::geom_point(size = 1.7) +
    ggplot2::scale_colour_manual(values = line_colours, name = NULL) +
    ggplot2::scale_linetype_manual(values = c("solid", "22", "dotted")[seq_along(line_levels)], name = NULL) +
    ggplot2::labs(x = "Sampling effort", y = "Functional change") + common_theme
}
