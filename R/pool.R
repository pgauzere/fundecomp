#' Summarise species, abundance, and trait information across a selected pool
#'
#' The default pool is retrospective and uses all dates in the selected scope.
#' Each date contributes one normalized SAD before date-level averaging, so
#' dates are not silently weighted by their number of samples. `sad` is returned
#' as a descriptive pool summary. It does not define a regional-pool baseline
#' for `decompose_fmetric_change()` or `decompose_fmetric_ts()`; those
#' decompositions use only the observed reference and focal endpoints.
#'
#' @param x An fundecomp_input object.
#' @param pool_species Optional admissible species vector.
#' @param baseline Pool reference definition.
#' @param sad_pool Definition used only for the descriptive `sad` summary in
#'   the returned pool object. It does not alter direct endpoint decomposition.
#' @param trait_pool Trait pool definition.
#' @param pool_scope Pool scope; "dataset", "series", or "custom".
#' @param series Series used for pool_scope = "series".
#' @param time_reference A list with series and time used by
#'   baseline = "time_reference".
#' @param exclude_time Optional time value omitted when estimating the pooled
#'   SAD. Traits and the admissible species pool remain from the selected scope.
#' @param exclude_series Optional series identifier used with exclude_time.
#' @param custom_pool Optional list with species, sad, traits, and
#'   optional trait_weights.
#' @return An object of class fundecomp_pool.
#' @details The returned object is a descriptive summary. Current decomposition
#'   functions do not use it as a counterfactual baseline; a non-`NULL` `pool`
#'   supplied to a decomposition is rejected.
#' @export
make_fundecomp_pool <- function(x, pool_species = NULL,
                           baseline = c("dataset_pool", "time_reference", "custom_pool"),
                           sad_pool = c("pooled_empirical", "uniform", "custom"),
                           trait_pool = c("pooled_species", "abundance_weighted", "custom"),
                           pool_scope = c("dataset", "series", "custom"),
                           series = NULL, time_reference = NULL,
                           exclude_time = NULL, exclude_series = NULL,
                           custom_pool = NULL) {
  if (!inherits(x, "fundecomp_input")) .fd_abort("x must be an fundecomp_input object.")
  baseline <- match.arg(baseline)
  sad_pool <- match.arg(sad_pool)
  trait_pool <- match.arg(trait_pool)
  pool_scope <- match.arg(pool_scope)
  needs_custom <- baseline == "custom_pool" || sad_pool == "custom" ||
    trait_pool == "custom" || pool_scope == "custom"
  if (needs_custom && (is.null(custom_pool) || !is.list(custom_pool))) {
    .fd_abort("custom_pool must be supplied for custom pool definitions.")
  }
  if (baseline == "time_reference") {
    if (!is.list(time_reference) || !all(c("series", "time") %in% names(time_reference))) {
      .fd_abort("time_reference must be a list with series and time for baseline = time_reference.")
    }
  }
  if (pool_scope == "series") {
    if (length(series) != 1L || is.na(series)) .fd_abort("series must identify one series for pool_scope = series.")
    series <- as.character(series)
  }

  all_states <- .fd_states(x)
  states <- all_states
  if (pool_scope == "series") states <- Filter(function(z) identical(z$series, series), states)
  if (pool_scope == "series" && !length(states)) .fd_abort("Requested pool series is absent from the input.")
  if (baseline == "time_reference") {
    states <- list(.fd_state(x, time_reference$series, time_reference$time))
  }
  if (pool_scope == "custom") states <- list()
  scope_states <- states

  if (!is.null(exclude_time)) {
    if (length(exclude_time) != 1L || is.na(exclude_time)) {
      .fd_abort("exclude_time must be one non-missing time value.")
    }
    if (!is.null(exclude_series) && length(exclude_series) != 1L) {
      .fd_abort("exclude_series must be one series identifier.")
    }
    omit <- vapply(states, function(z) {
      same_time <- .fd_time_key(z$time) == .fd_time_key(exclude_time)
      same_series <- is.null(exclude_series) || identical(z$series, as.character(exclude_series))
      same_time && same_series
    }, logical(1))
    if (!any(omit)) .fd_abort("exclude_time is absent from the selected pool scope.")
    states <- states[!omit]
  }
  if (!length(states) && sad_pool != "custom") {
    .fd_abort("No states remain to estimate the selected pool.")
  }

  # Leave-one-time-out changes the SAD estimation set, not the admissible
  # species pool. This avoids silently removing species occurring only at the
  # omitted date.
  observed <- unique(unlist(lapply(scope_states, function(z) colnames(z$comm))))
  custom_species <- if (is.list(custom_pool)) custom_pool$species else NULL
  if (is.null(pool_species)) {
    species <- if (!is.null(custom_species)) as.character(custom_species) else observed
    if (!length(species) && is.list(custom_pool) && !is.null(custom_pool$traits)) {
      species <- rownames(as.data.frame(custom_pool$traits))
    }
  } else species <- as.character(pool_species)
  if (!length(species) || anyNA(species) || any(!nzchar(species)) || anyDuplicated(species)) {
    .fd_abort("pool_species must be a non-empty unique species vector.")
  }

  if (trait_pool == "custom") {
    trait_table <- as.data.frame(custom_pool$traits, stringsAsFactors = FALSE)
    if (is.null(rownames(trait_table))) .fd_abort("Custom traits must have species row names.")
    if (!all(species %in% rownames(trait_table))) .fd_abort("Custom traits lacks selected pool species.")
    trait_table <- trait_table[species, , drop = FALSE]
    if (!ncol(trait_table) || !all(vapply(trait_table, is.numeric, logical(1))) ||
        any(!is.finite(as.matrix(trait_table)))) {
      .fd_abort("Custom pool traits must be finite numeric columns.")
    }
  } else {
    if (!all(species %in% rownames(x$traits))) .fd_abort("The selected pool contains species without complete traits.")
    trait_table <- x$traits[species, , drop = FALSE]
  }
  rownames(trait_table) <- species

  state_sad <- lapply(states, function(z) {
    a <- colSums(z$comm)
    out <- stats::setNames(numeric(length(species)), species)
    common <- intersect(names(a), species)
    out[common] <- a[common]
    if (sum(out) > 0) out / sum(out) else out
  })
  if (sad_pool == "pooled_empirical") {
    sad <- Reduce("+", lapply(state_sad, function(p) p / length(state_sad)))
    if (sum(sad) <= 0) .fd_abort("The selected pool has no positive abundance.")
    sad <- sad / sum(sad)
  } else if (sad_pool == "uniform") {
    sad <- stats::setNames(rep(1 / length(species), length(species)), species)
  } else {
    sad <- custom_pool$sad
    if (is.null(names(sad))) .fd_abort("Custom sad must be a named numeric vector.")
    sad <- as.numeric(sad[species])
    if (anyNA(sad) || any(sad < 0) || sum(sad) <= 0) .fd_abort("Custom sad is invalid for the selected pool.")
    sad <- stats::setNames(sad / sum(sad), species)
  }

  if (trait_pool == "abundance_weighted") {
    weights <- sad
  } else if (trait_pool == "pooled_species") {
    weights <- stats::setNames(rep(1 / length(species), length(species)), species)
  } else if (!is.null(custom_pool$trait_weights)) {
    weights <- as.numeric(custom_pool$trait_weights[species])
    if (anyNA(weights) || any(weights < 0) || sum(weights) <= 0) .fd_abort("Custom trait_weights is invalid.")
    weights <- stats::setNames(weights / sum(weights), species)
  } else {
    weights <- stats::setNames(rep(1 / length(species), length(species)), species)
  }
  attr(trait_table, "weights") <- weights
  baseline_name <- baseline
  if (!is.null(exclude_time)) baseline_name <- paste0(baseline_name, "_leave_one_time_out")
  structure(list(
    species = species,
    sad = sad,
    traits = trait_table,
    trait_names = names(trait_table),
    trait_weights = weights,
    baseline = baseline_name,
    sad_pool = sad_pool,
    trait_pool = trait_pool,
    pool_scope = pool_scope,
    states = states,
    diagnostics = list(
      n_states = length(states), retrospective = baseline == "dataset_pool",
      excluded_time = exclude_time, excluded_series = exclude_series,
      reference = time_reference, series = series,
      sad_role = "descriptive_pool_summary_only",
      direct_decomposition_sad = "observed_reference_and_focal_endpoints"
    )
  ), class = "fundecomp_pool")
}
