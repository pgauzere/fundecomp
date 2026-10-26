#' Decompose a metric across a temporal series
#'
#' @param x An `fundecomp_input` object.
#' @param metric Metric name or custom function.
#' @param series Optional series identifiers.
#' @param time_reference `"first"` (the validated default), `"previous"`, or
#'   one time value per selected series.
#' @param focal_times Optional focal times; defaults to all dates after the
#'   reference date.
#' @param ... Arguments passed to `decompose_fmetric_change()` for each pair.
#' @return An object of class `fundecomp_timeseries` containing a long table and per-change
#'   `fundecomp_result` objects.
#' @export
decompose_fmetric_ts <- function(x, metric = "FRic", series = NULL,
                                 time_reference = c("first", "previous"),
                                 focal_times = NULL, ...) {
  if (!inherits(x, "fundecomp_input")) .fd_abort("`x` must be an `fundecomp_input` object.")
  if (!isTRUE(x$has_time)) .fd_abort("`decompose_fmetric_ts()` requires an explicitly supplied ordered `time` variable.")
  dots <- list(...)
  if (length(time_reference) == 2L && identical(time_reference, c("first", "previous"))) {
    time_reference <- "first"
  } else if (length(time_reference) == 1L && time_reference %in% c("first", "previous")) {
    time_reference <- time_reference
  } else if (length(time_reference) != 1L) {
    .fd_abort("`time_reference` must be `\"first\"`, `\"previous\"`, or one supplied time value.")
  }
  all_series <- unique(x$meta$series)
  if (is.null(series)) series <- all_series
  if (!all(series %in% all_series)) .fd_abort("Unknown series in `series`.")
  outputs <- list()
  counter <- 0L
  for (s in series) {
    idx <- x$meta$series == s
    keys <- unique(.fd_time_key(x$meta$time[idx]))
    first_idx <- match(keys, .fd_time_key(x$meta$time[idx]))
    times <- x$meta$time[idx][first_idx]
    times <- times[.fd_order_time_values(times)]
    if (length(times) < 2L) next
    if (identical(time_reference, "first")) ref_i <- 1L else if (identical(time_reference, "previous")) ref_i <- NA_integer_ else {
      ref_i <- match(.fd_time_key(time_reference), .fd_time_key(times))
      if (is.na(ref_i)) .fd_abort("Requested time reference is absent from series `", s, "`.")
    }
    chosen <- if (is.null(focal_times)) seq_along(times) else {
      match(.fd_time_key(focal_times), .fd_time_key(times))
    }
    if (anyNA(chosen)) .fd_abort("One or more `focal_times` are absent from series `", s, "`.")
    for (j in chosen) {
      ri <- if (is.na(ref_i)) j - 1L else ref_i
      if (ri < 1L || ri == j) next
      counter <- counter + 1L
      args <- c(list(x = x,
                     reference = list(series = s, time = times[ri]),
                     focal = list(series = s, time = times[j]),
                     metric = metric), dots)
      outputs[[counter]] <- do.call(decompose_fmetric_change, args)
    }
  }
  if (!length(outputs)) .fd_abort("No valid temporal comparisons were selected.")
  curves <- do.call(rbind, lapply(outputs, function(z) z$curves))
  .fd_time_delta <- function(a, b) {
    if (inherits(a, "POSIXt") || inherits(b, "POSIXt")) {
      return(as.numeric(difftime(b, a, units = "secs")))
    }
    if (inherits(a, "Date") || inherits(b, "Date")) return(as.numeric(b - a))
    if (is.numeric(a) && is.numeric(b)) return(as.numeric(b) - as.numeric(a))
    NA_real_
  }
  deltas <- vapply(outputs, function(z) .fd_time_delta(z$reference$time, z$focal$time), numeric(1))
  delta_units <- if (inherits(x$meta$time, "POSIXt")) "seconds" else if (inherits(x$meta$time, "Date")) "days" else if (is.numeric(x$meta$time)) "time_units" else NA_character_
  curves$delta_t <- rep(deltas, times = vapply(outputs, function(z) nrow(z$curves), integer(1)))
  curves$delta_t_unit <- rep(delta_units, nrow(curves))
  out <- list(curves = curves, changes = outputs,
              time_reference = time_reference, metric = metric)
  class(out) <- "fundecomp_timeseries"
  out
}
