.fd_metric_registry <- new.env(parent = emptyenv())

#' Register a functional metric
#'
#' @param name Metric name.
#' @param fun Function accepting `abundance`, `traits`, and `trait`.
#' @param required_traits Minimum number of trait dimensions.
#' @param effort_unit Effort scales supported by the metric.
#' @param min_species Minimum number of observed species required.
#' @export
register_fmetric <- function(name, fun, required_traits = 1L,
                             effort_unit = c("individual", "sample", "full"),
                             min_species = 1L) {
  if (length(name) != 1L || !nzchar(name) || !is.function(fun)) {
    .fd_abort("`name` must be non-empty and `fun` must be a function.")
  }
  .fd_metric_registry[[name]] <- list(
    name = name, fun = fun, required_traits = as.integer(required_traits),
    effort_unit = effort_unit, min_species = as.integer(min_species)
  )
  invisible(name)
}

#' List registered functional metrics
#' @export
available_fmetrics <- function() sort(ls(.fd_metric_registry))

.fd_get_metric <- function(metric) {
  legacy <- c(biomass = "trait_sum", variance = "CWV")
  if (is.character(metric) && length(metric) == 1L && metric %in% names(legacy)) {
    metric <- unname(legacy[metric])
  }
  if (is.list(metric) && is.function(metric$fun)) return(metric)
  if (is.function(metric)) return(list(name = "custom", fun = metric, required_traits = 1L,
                                       effort_unit = c("individual", "sample", "full"), min_species = 1L))
  if (length(metric) != 1L || !exists(metric, .fd_metric_registry, inherits = FALSE)) {
    .fd_abort("Unknown metric `", metric, "`. Available metrics: ", paste(available_fmetrics(), collapse = ", "), ".")
  }
  .fd_metric_registry[[metric]]
}

.fd_trait_matrix <- function(traits, trait = NULL, required = 1L) {
  if (is.null(traits) || !nrow(traits)) .fd_abort("No trait data are available.")
  trait_names <- setdiff(names(traits), attr(traits, "species_col"))
  if (!length(trait_names)) trait_names <- names(traits)
  if (is.null(trait)) trait <- trait_names
  if (is.numeric(trait)) trait <- trait_names[trait]
  if (!all(trait %in% trait_names)) .fd_abort("Requested trait is not in the trait table.")
  if (length(trait) < required) .fd_abort("Metric requires at least ", required, " trait dimension(s).")
  out <- as.matrix(traits[, trait, drop = FALSE])
  storage.mode(out) <- "double"
  if (any(!is.finite(out))) .fd_abort("Traits must be finite for metric calculation.")
  out
}

.fd_metric_biomass <- function(abundance, traits, trait = NULL, ...) {
  tr <- .fd_single_trait(traits, trait)
  sum(abundance * tr)
}

.fd_metric_cwm <- function(abundance, traits, trait = NULL, ...) {
  tr <- .fd_single_trait(traits, trait)
  if (sum(abundance) <= 0) return(NA_real_)
  sum(abundance * tr) / sum(abundance)
}

.fd_metric_cwv <- function(abundance, traits, trait = NULL, ...) {
  tr <- .fd_single_trait(traits, trait)
  if (sum(abundance) <= 0) return(NA_real_)
  mu <- sum(abundance * tr) / sum(abundance)
  sum(abundance * (tr - mu)^2) / sum(abundance)
}

.fd_single_trait <- function(traits, trait = NULL) {
  available <- setdiff(names(traits), attr(traits, "species_col"))
  if (is.null(trait)) {
    if (length(available) != 1L) {
      .fd_abort("This metric requires one trait; supply `trait` when multiple traits are available.")
    }
    trait <- available
  }
  if (length(trait) != 1L) .fd_abort("This metric requires exactly one trait.")
  .fd_trait_matrix(traits, trait, 1L)[, 1L]
}

.fd_metric_fric <- function(abundance, traits, trait = NULL, ...) {
  tr <- .fd_trait_matrix(traits, trait, 1L)
  present <- abundance > 0
  tr <- tr[present, , drop = FALSE]
  if (nrow(tr) < 2L) return(0)
  n_dimensions <- ncol(tr)
  if (n_dimensions == 1L) return(diff(range(tr[, 1L])))

  # Multivariate FRic is the volume of the convex hull in the supplied trait
  # coordinate system. It is deliberately unweighted by abundance. A hull with
  # fewer than d + 1 points, or one contained in a lower-dimensional affine
  # subspace, has zero d-dimensional volume.
  centred <- sweep(tr, 2L, colMeans(tr), "-")
  if (nrow(tr) <= n_dimensions || qr(centred, tol = sqrt(.Machine$double.eps))$rank < n_dimensions) {
    return(0)
  }
  if (!requireNamespace("geometry", quietly = TRUE)) {
    .fd_abort(
      "Multivariate FRic requires the suggested package `geometry`. ",
      "Install it with install.packages(\"geometry\")."
    )
  }
  hull <- tryCatch(
    geometry::convhulln(tr, options = "Qt", output.options = "FA"),
    error = function(e) NULL
  )
  volume <- if (is.null(hull)) NA_real_ else as.numeric(hull$vol)[1L]
  if (!is.finite(volume)) return(0)
  max(0, volume)
}

.fd_metric_fdis <- function(abundance, traits, trait = NULL, ...) {
  tr <- .fd_trait_matrix(traits, trait, 1L)
  present <- abundance > 0
  if (!any(present) || sum(abundance[present]) <= 0) return(NA_real_)
  tr <- tr[present, , drop = FALSE]
  a <- abundance[present]
  centroid <- colSums(tr * a) / sum(a)
  d <- sqrt(rowSums((tr - matrix(centroid, nrow(tr), ncol(tr), byrow = TRUE))^2))
  sum(a * d) / sum(a)
}

.fd_metric_richness <- function(abundance, traits = NULL, ...) sum(abundance > 0)

.fd_register_defaults <- function() {
  register_fmetric("richness", .fd_metric_richness, 0L,
                   c("individual", "sample", "full"), 1L)
  register_fmetric("trait_sum", .fd_metric_biomass, 1L, c("individual", "sample", "full"), 1L)
  # Compatibility aliases remain discoverable for v0.1 callers. The resolver
  # warns and redirects them to the names used by current analyses.
  register_fmetric("biomass", .fd_metric_biomass, 1L, c("individual", "sample", "full"), 1L)
  register_fmetric("CWM", .fd_metric_cwm, 1L, c("individual", "sample", "full"), 1L)
  register_fmetric("CWV", .fd_metric_cwv, 1L, c("individual", "sample", "full"), 1L)
  register_fmetric("variance", .fd_metric_cwv, 1L, c("individual", "sample", "full"), 1L)
  register_fmetric("FRic", .fd_metric_fric, 1L, c("individual", "sample", "full"), 1L)
  register_fmetric("FDis", .fd_metric_fdis, 1L, c("individual", "sample", "full"), 1L)
}
.fd_register_defaults()

.fd_sample_counts <- function(abundance, m) {
  nms <- names(abundance)
  abundance <- as.numeric(abundance)
  names(abundance) <- nms
  N <- sum(abundance)
  if (m < 0 || m > N) .fd_abort("Individual effort must be between zero and total abundance.")
  if (m == 0) return(stats::setNames(rep(0, length(abundance)), names(abundance)))
  if (m == N) return(stats::setNames(abundance, names(abundance)))
  out <- numeric(length(abundance))
  remaining_m <- as.numeric(m)
  remaining_N <- as.numeric(N)
  if (length(abundance) > 1L) for (i in seq_len(length(abundance) - 1L)) {
    take <- if (remaining_m <= 0 || abundance[i] <= 0) 0 else
      stats::rhyper(1L, abundance[i], remaining_N - abundance[i], remaining_m)
    out[i] <- take
    remaining_m <- remaining_m - take
    remaining_N <- remaining_N - abundance[i]
  }
  out[length(abundance)] <- remaining_m
  stats::setNames(out, names(abundance))
}

.fd_aggregate_rows <- function(comm, rows) {
  a <- colSums(comm[rows, , drop = FALSE])
  stats::setNames(as.numeric(a), colnames(comm))
}

.fd_spatial_orders <- function(state, spat_algo = "kNN") {
  if (is.null(state$coord_names) || length(state$coord_names) != 2L) {
    .fd_abort("Spatial sample-based effort requires two coordinate columns.")
  }
  coords <- as.matrix(state$meta[state$coord_names])
  if (any(!is.finite(coords))) .fd_abort("Spatial coordinates must be finite.")
  spat_algo <- .fd_check_choice(spat_algo, c("kNN", "kNCN"), "spat_algo")
  if (spat_algo == "kNCN" && isTRUE(state$latlong)) {
    .fd_abort(
      '`spat_algo = "kNCN"` is experimental and is not supported with ',
      "longitude/latitude coordinates. Use `spat_algo = \"kNN\"` or project ",
      "coordinates to a suitable planar coordinate system."
    )
  }
  n <- nrow(coords)
  dmat <- if (isTRUE(state$latlong)) .fd_sphere_dist(coords) else as.matrix(stats::dist(coords))
  out <- vector("list", n)
  for (focal in seq_len(n)) {
    if (spat_algo == "kNN") {
      out[[focal]] <- order(dmat[focal, ], seq_len(n))
    } else {
      chosen <- focal
      remaining <- setdiff(seq_len(n), chosen)
      while (length(remaining)) {
        centroid <- colMeans(coords[chosen, , drop = FALSE])
        d <- rowSums((coords[remaining, , drop = FALSE] -
                        matrix(centroid, nrow = length(remaining), ncol = 2L, byrow = TRUE))^2)
        next_plot <- remaining[order(d, remaining)[1L]]
        chosen <- c(chosen, next_plot)
        remaining <- setdiff(remaining, next_plot)
      }
      out[[focal]] <- chosen
    }
  }
  out
}

.fd_sphere_dist <- function(coords, r = 6378.137) {
  # Same Haversine implementation as mobr::sphere_dist, kept local because
  # sphere_dist is internal to mobr. Coordinates are longitude, latitude.
  coslat1 <- cos((coords[, 2] * pi) / 180)
  sinlat1 <- sin((coords[, 2] * pi) / 180)
  coslon1 <- cos((coords[, 1] * pi) / 180)
  sinlon1 <- sin((coords[, 1] * pi) / 180)
  xyz <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1)
  pp <- xyz %*% t(xyz)
  r * acos(ifelse(abs(pp) > 1, sign(pp), pp))
}

.fd_expand_itv_support <- function(abundance, traits, itv) {
  values <- list()
  weights <- list()
  labels <- character()
  for (sp in names(abundance)) {
    residuals <- itv[[sp]]
    if (is.null(residuals) || !nrow(residuals)) {
      residuals <- matrix(0, nrow = 1L, ncol = ncol(traits),
                          dimnames = list(NULL, names(traits)))
    }
    residuals <- as.matrix(residuals[, names(traits), drop = FALSE])
    values[[sp]] <- sweep(residuals, 2L, as.numeric(traits[sp, ]), "+")
    weights[[sp]] <- rep(as.numeric(abundance[sp]) / nrow(residuals), nrow(residuals))
    labels <- c(labels, paste0(sp, "::trait_observation_", seq_len(nrow(residuals))))
  }
  out_traits <- as.data.frame(do.call(rbind, values), stringsAsFactors = FALSE)
  rownames(out_traits) <- labels
  list(abundance = stats::setNames(unlist(weights, use.names = FALSE), labels),
       traits = out_traits)
}

.fd_eval_abundance <- function(abundance, traits, metric, trait = NULL, itv = NULL) {
  keep <- abundance > 0
  abundance <- abundance[keep]
  traits <- traits[names(abundance), , drop = FALSE]
  if (!length(abundance)) return(NA_real_)
  if (identical(metric$name, "richness")) return(sum(abundance > 0))
  if (!is.null(itv)) {
    support <- .fd_expand_itv_support(abundance, traits, itv[names(abundance)])
    abundance <- support$abundance
    traits <- support$traits
  }
  if (length(abundance) < metric$min_species) return(NA_real_)
  metric$fun(abundance = abundance, traits = traits, trait = trait)
}

.fd_expected_richness_ibr <- function(abundance, effort) {
  abundance <- as.numeric(abundance)
  N <- sum(abundance)
  if (effort == 0 || N == 0) return(0)
  if (effort == N) return(sum(abundance > 0))
  present <- abundance > 0
  a <- abundance[present]
  log_denom <- lchoose(N, effort)
  ratios <- vapply(a, function(ai) {
    if (N - ai < effort) return(0)
    exp(lchoose(N - ai, effort) - log_denom)
  }, numeric(1))
  sum(-expm1(log(ratios)), na.rm = TRUE)
}

.fd_eval_curve_state <- function(state, metric, effort = NULL, effort_unit = "full",
                                   trait = NULL, n_rep = 1L, seed = NULL,
                                   spat_algo = NULL, return_replicates = FALSE) {
  .fd_with_seed(seed, function() .fd_eval_curve_state_impl(
    state = state, metric = metric, effort = effort, effort_unit = effort_unit,
    trait = trait, n_rep = n_rep, spat_algo = spat_algo,
    return_replicates = return_replicates))
}

.fd_eval_curve_state_impl <- function(state, metric, effort = NULL,
                                        effort_unit = "full", trait = NULL,
                                        n_rep = 1L, spat_algo = NULL,
                                        return_replicates = FALSE) {
  metric <- .fd_get_metric(metric)
  effort_unit <- .fd_check_choice(effort_unit, c("individual", "sample", "full"), "effort_unit")
  if (!effort_unit %in% metric$effort_unit) .fd_abort("Metric does not support effort unit `", effort_unit, "`.")
  if (length(n_rep) != 1L || !is.finite(n_rep) || n_rep < 1L || n_rep != floor(n_rep)) {
    .fd_abort("`n_rep` must be one positive integer.")
  }
  n_rep <- as.integer(n_rep)
  comm <- state$comm
  traits <- state$traits
  itv <- state$itv
  if (is.null(names(traits))) rownames(traits) <- colnames(comm)
  traits <- traits[colnames(comm), , drop = FALSE]
  if (effort_unit == "full") {
    value <- .fd_eval_abundance(colSums(comm), traits, metric, trait, itv)
    reps <- matrix(rep(value, n_rep), nrow = n_rep, ncol = 1L)
    colnames(reps) <- "0"
    return(if (return_replicates) reps else value)
  }
  if (is.null(effort)) .fd_abort("`effort` is required unless `effort_unit = \"full\"`.")
  if (!length(effort) || any(!is.finite(effort)) || any(effort < 0)) {
    .fd_abort("`effort` must contain non-negative finite values.")
  }
  if (any(effort != floor(effort))) .fd_abort("Individual and sample effort values must be integers.")
  effort <- as.integer(effort)
  sorted_effort <- sort(unique(effort))
  effort_map <- match(effort, sorted_effort)
  if (effort_unit == "individual") {
    if (!isTRUE(state$count_data) || any(comm != floor(comm))) {
      .fd_abort("Individual rarefaction requires integer count data; use sample effort for continuous abundances.")
    }
    if (any(sorted_effort > sum(comm))) .fd_abort("Individual effort exceeds total abundance.")
    abundance <- colSums(comm)
    if (identical(metric$name, "richness")) {
      values <- vapply(sorted_effort, function(m)
        .fd_expected_richness_ibr(abundance, m), numeric(1))
      reps <- matrix(rep(values, each = n_rep), nrow = n_rep)
    } else {
      reps <- matrix(NA_real_, nrow = n_rep, ncol = length(sorted_effort))
      for (r in seq_len(n_rep)) {
        current <- stats::setNames(rep(0, length(abundance)), names(abundance))
        remaining <- abundance
        previous <- 0L
        for (j in seq_along(sorted_effort)) {
          increment <- sorted_effort[j] - previous
          if (increment > 0L) {
            draw <- .fd_sample_counts(remaining, increment)
            current <- current + draw
            remaining <- remaining - draw
          }
          reps[r, j] <- .fd_eval_abundance(current, traits, metric, trait, itv)
          previous <- sorted_effort[j]
        }
      }
    }
  } else {
    if (any(sorted_effort > nrow(comm))) .fd_abort("Sample effort exceeds the number of samples.")
    if (!is.null(spat_algo)) {
      orders <- .fd_spatial_orders(state, spat_algo)
      spatial_reps <- matrix(NA_real_, nrow = length(orders), ncol = length(sorted_effort))
      for (r in seq_along(orders)) {
        ord <- orders[[r]]
        for (j in seq_along(sorted_effort)) {
          rows <- if (sorted_effort[j] == 0L) integer() else ord[seq_len(sorted_effort[j])]
          a <- if (length(rows)) .fd_aggregate_rows(comm, rows) else
            stats::setNames(rep(0, ncol(comm)), colnames(comm))
          spatial_reps[r, j] <- .fd_eval_abundance(a, traits, metric, trait, itv)
        }
      }
      values <- colMeans(spatial_reps)
      reps <- matrix(rep(values, each = n_rep), nrow = n_rep)
    } else {
      reps <- matrix(NA_real_, nrow = n_rep, ncol = length(sorted_effort))
      for (r in seq_len(n_rep)) {
        ord <- sample.int(nrow(comm), nrow(comm), replace = FALSE)
        for (j in seq_along(sorted_effort)) {
          rows <- if (sorted_effort[j] == 0L) integer() else ord[seq_len(sorted_effort[j])]
          a <- if (length(rows)) .fd_aggregate_rows(comm, rows) else
            stats::setNames(rep(0, ncol(comm)), colnames(comm))
          reps[r, j] <- .fd_eval_abundance(a, traits, metric, trait, itv)
        }
      }
    }
  }
  reps <- reps[, effort_map, drop = FALSE]
  colnames(reps) <- as.character(effort)
  if (return_replicates) reps else colMeans(reps)
}

.fd_eval_state <- function(state, metric, effort = NULL, effort_unit = "full",
                             trait = NULL, n_rep = 1L, seed = NULL,
                             spat_algo = NULL) {
  value <- .fd_eval_curve_state(state, metric, effort, effort_unit, trait,
                                  n_rep, seed, spat_algo, FALSE)
  if (length(value) != 1L) .fd_abort("Internal state evaluation expected one effort value.")
  as.numeric(value)
}

#' Evaluate a metric for one observed state over an effort grid
#' @param x An `fundecomp_input` object or explicit internal state list.
#' @param state State identifier when `x` contains multiple states. Omission is
#'   allowed only when exactly one state is present.
#' @param metric Metric name or function.
#' @param effort Effort grid.
#' @param effort_unit `"individual"`, `"sample"`, or `"full"`.
#' @param trait Trait column(s) used by the metric.
#' @param n_rep Number of Monte Carlo replicates. Within each replicate, effort
#'   points are nested prefixes of one without-replacement individual draw or
#'   one sample order.
#' @param seed Random seed. When `NULL`, a local effective seed is generated,
#'   returned in the output, and does not alter the caller's RNG state.
#' @param spat_algo Optional spatial ordering. `"kNN"` is the supported
#'   algorithm. `"kNCN"` is experimental and is rejected for longitude/latitude
#'   coordinates.
#' @param series Optional series identifier when state labels repeat across series.
#' @return A data frame with effort and metric value. The matched replicate
#'   matrix is stored in the `replicates` attribute and `nested_sampling`
#'   records whether effort points share a nested draw.
#' @export
rarefy_community <- function(x, state = NULL, metric = "FRic", effort = NULL, effort_unit = "individual",
                             trait = NULL, n_rep = 1L, seed = NULL,
                             spat_algo = NULL, series = NULL) {
  if (inherits(x, "fundecomp_input")) {
    if (length(metric) > 1L && is.character(state) && length(state) == 1L) {
      old_metric <- state; old_effort <- metric; old_unit <- effort
      state <- NULL; metric <- old_metric; effort <- old_effort
      if (is.character(old_unit) && length(old_unit) == 1L) effort_unit <- old_unit
    }
    all_states <- .fd_states(x)
    if (!is.null(series)) all_states <- Filter(function(z) z$series %in% as.character(series), all_states)
    if (!length(all_states)) .fd_abort("No states match the requested `series`.")
    if (is.null(state)) {
      if (length(all_states) != 1L) .fd_abort("`state` must be specified when `x` contains multiple states.")
      st <- all_states[[1L]]
    } else if (is.list(state) && !is.null(state$comm)) {
      st <- state
    } else {
      hit <- which(vapply(all_states, function(z) .fd_time_key(z$time) == .fd_time_key(state), logical(1)))
      if (length(hit) != 1L) .fd_abort("`state` must identify exactly one state in `x`.")
      st <- all_states[[hit]]
    }
  } else if (is.list(x) && !is.null(x$comm)) {
    if (is.function(state)) {
      metric <- state; state <- NULL
    } else if (is.character(state) && length(state) == 1L) {
      old_metric <- state
      if (!identical(metric, "FRic")) {
        old_effort <- metric; old_unit <- effort
        metric <- old_metric; effort <- old_effort
        if (is.character(old_unit) && length(old_unit) == 1L) effort_unit <- old_unit
      } else metric <- old_metric
    }
    st <- x
  } else .fd_abort("`x` must be an `fundecomp_input` object or an explicit state list.")
  metric <- .fd_get_metric(metric)
  if (effort_unit == "full") effort <- 0
  if (!length(effort) || any(!is.finite(effort)) || any(effort < 0)) .fd_abort("`effort` must be finite and non-negative.")
  if (effort_unit == "individual" && any(effort != floor(effort))) {
    .fd_abort("Individual effort must be integer.")
  }
  if (length(n_rep) != 1L || !is.finite(n_rep) || n_rep < 1L || n_rep != floor(n_rep)) {
    .fd_abort("`n_rep` must be a positive integer.")
  }
  effective_seed <- .fd_effective_seed(seed)
  reps <- .fd_eval_curve_state(st, metric, effort, effort_unit, trait,
                                 as.integer(n_rep), effective_seed, spat_algo, TRUE)
  vals <- colMeans(reps)
  seed_out <- effective_seed
  out <- data.frame(effort = effort, metric = metric$name, value = vals,
                    effort_unit = effort_unit, n_rep = as.integer(n_rep), seed = seed_out,
                    stringsAsFactors = FALSE)
  attr(out, "replicates") <- reps
  attr(out, "nested_sampling") <- effort_unit %in% c("individual", "sample")
  out
}
