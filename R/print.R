#' @export
print.fundecomp_input <- function(x, ...) {
  cat("fundecomp input\n")
  cat("  samples:", nrow(x$comm), " species:", ncol(x$comm),
      " series:", length(unique(x$meta$series)), "\n")
  cat("  missing traits dropped:", length(x$diagnostics$dropped_species), "\n")
  if (isTRUE(x$has_itv)) {
    cat("  trait observations:", x$diagnostics$n_trait_observations,
        " date-specific:", isTRUE(x$diagnostics$trait_observations_scoped), "\n")
  }
  invisible(x)
}

#' @export
print.fundecomp_pool <- function(x, ...) {
  cat("fundecomp pool\n")
  cat("  species:", length(x$species), " SAD:", x$sad_pool,
      " traits:", x$trait_pool, " scope:", x$pool_scope, "\n")
  invisible(x)
}

#' @export
print.fundecomp_result <- function(x, ...) {
  cat("fundecomp decomposition\n")
  cat("  metric:", x$metric, " method:", x$method,
      " effort unit:", x$effort_unit, "\n")
  cat("  reference:", x$reference$series, x$reference$time,
      " focal:", x$focal$time, "\n")
  if (!isTRUE(x$full_decomposition)) {
    cat("  selected components:", paste(unname(.fd_component_labels[x$components]),
      collapse = ", "), "\n")
    cat("  mean partial change:", format(mean(x$curves$partial_change), digits = 6),
        " observed total:", format(mean(x$curves$observed_total_change), digits = 6),
        " remaining difference:", format(mean(x$curves$remaining_difference), digits = 6), "\n")
  }
  if (!is.null(x$ta_axis)) {
    axis_label <- if (identical(x$ta_axis$type, "trait_column")) {
      paste0("trait column `", x$ta_axis$name, "`")
    } else "named species scores"
    cat("  T-A axis:", axis_label, "\n")
  }
  cat("  max reconstruction error:", format(max(abs(x$curves$reconstruction_error)), scientific = TRUE), "\n")
  if ("component_label" %in% names(x$curves)) {
    cat("  components:", paste(unique(x$curves$component_label), collapse = ", "), "\n")
  }
  invisible(x)
}

#' @export
summary.fundecomp_result <- function(object, ...) {
  list(metric = object$metric, method = object$method,
       reference = object$reference, focal = object$focal,
       components = object$components, full_decomposition = object$full_decomposition,
       partial_change = unique(object$curves$partial_change),
       observed_total_change = unique(object$curves$observed_total_change),
       remaining_difference = unique(object$curves$remaining_difference),
       ta_axis = object$ta_axis,
       stats = get_fdelta_stats(object),
       max_reconstruction_error = max(abs(object$curves$reconstruction_error)),
       constraints = object$constraints)
}
