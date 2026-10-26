#' Construct a validated temporal functional community object
#'
#' @param comm Community matrix, with samples in rows and species in columns.
#' @param long Optional long table with `series`, `state` (or legacy `time`), `sample`, `species`,
#'   and `abundance` columns.
#' @param traits Species-by-trait table.
#' @param trait_observations Optional long table of individual trait
#'   measurements. It must contain `species` and every trait column; supplying
#'   both `series` and `time` makes the within-species distributions date
#'   specific.
#' @param plot_attr Optional sample metadata, including coordinates.
#' @param series,state,sample Optional vectors identifying generic assemblage states.
#' @param time Optional ordered temporal variable; use only for temporal analyses.
#' @param species_col Species identifier column in `traits`.
#' @param trait_na Missing-trait policy: `"drop"` or `"error"`; dropped
#'   species are returned in diagnostics.
#' @param itv_na Policy for a species present at a date but lacking a
#'   date-specific trait observation: `"error"` (default) or
#'   `"species_mean"`, which uses the species-level trait value with zero
#'   within-species deviation.
#' @param coord_names Optional names of the two coordinate columns.
#' @param latlong Whether the coordinates are longitude/latitude.
#' @param count_data Whether abundances must be integer counts.
#' @param exposure Optional exposure or effort vector, retained as metadata.
#' @param options Additional options recorded in the object.
#' @return An object of class `fundecomp_input`.
#' @export
make_fundecomp_input <- function(comm = NULL, long = NULL, traits, plot_attr = NULL,
                         series = NULL, time = NULL, sample = NULL, state = NULL,
                         species_col = NULL, trait_na = c("drop", "error"),
                         trait_observations = NULL,
                         itv_na = c("error", "species_mean"),
                         coord_names = NULL, latlong = FALSE,
                         count_data = TRUE, exposure = NULL, options = list()) {
  trait_na <- match.arg(trait_na)
  itv_na <- match.arg(itv_na)
  has_time <- !is.null(time)
  if (!is.null(comm) && !is.null(long)) .fd_abort("Supply either `comm` or `long`, not both.")
  if (is.null(comm) && is.null(long)) .fd_abort("Supply a community matrix or a long community table.")

  if (!is.null(long)) {
    long <- as.data.frame(long, stringsAsFactors = FALSE)
    if (!"state" %in% names(long) && "time" %in% names(long)) {
      long$state <- long$time
      has_time <- TRUE
    }
    required <- c("series", "state", "sample", "species", "abundance")
    if (!all(required %in% names(long))) {
      .fd_abort("Long input must contain: ", paste(required, collapse = ", "), ".")
    }
    if (anyDuplicated(long[c("series", "state", "sample", "species")])) {
      .fd_abort("Long input must have one row per series-time-sample-species combination.")
    }
    if (anyNA(long$species) || any(!nzchar(as.character(long$species)))) {
      .fd_abort("Long input species identifiers must be non-missing and non-empty.")
    }
    species <- unique(as.character(long$species))
    samples <- unique(long[c("series", "state", "sample")])
    sample_keys <- paste(samples$series, samples$state, samples$sample, sep = "@@")
    comm <- matrix(0, nrow = nrow(samples), ncol = length(species),
                   dimnames = list(sample_keys, species))
    row_id <- match(paste(long$series, long$state, long$sample, sep = "@@"),
                    paste(samples$series, samples$state, samples$sample, sep = "@@"))
    col_id <- match(as.character(long$species), species)
    comm[cbind(row_id, col_id)] <- long$abundance
    series <- samples$series
    state <- samples$state
    sample <- samples$sample
    plot_attr <- samples
    plot_attr <- plot_attr[c("series", "state", "sample")]
    extra <- setdiff(names(long), required)
    if (length(extra)) {
      extra_rows <- unique(long[c("series", "state", "sample", extra)])
      if (anyDuplicated(extra_rows[c("series", "state", "sample")])) {
        .fd_abort("Long metadata must be constant within each series-time-sample combination.")
      }
      extra_rows <- extra_rows[match(paste(samples$series, samples$state, samples$sample, sep = "@@"),
                                     paste(extra_rows$series, extra_rows$state, extra_rows$sample, sep = "@@")), , drop = FALSE]
      plot_attr <- cbind(plot_attr, extra_rows[extra])
    }
  }

  comm <- .fd_as_matrix(comm)
  .fd_check_abundance(comm, count_data)
  if (!is.null(state)) {
    if (!is.null(time)) .fd_abort("Supply `state` or `time`, not both.")
    time <- state
  }
  if (is.null(time) && !is.null(plot_attr) && "state" %in% names(plot_attr)) time <- plot_attr$state
  meta <- .fd_make_meta(nrow(comm), series, time, sample, plot_attr)
  meta$state <- meta$time
  if (!is.null(exposure)) {
    if (length(exposure) != nrow(comm) || any(!is.finite(exposure)) || any(exposure <= 0)) {
      .fd_abort("`exposure` must contain one positive finite value per community row.")
    }
    meta$exposure <- as.numeric(exposure)
  }
  coords <- .fd_coordinate_names(meta, coord_names)
  tt <- .fd_trait_table(traits, species_col)
  trait_data <- tt$data
  observed_species <- colnames(comm)[colSums(comm) > 0]
  absent_traits <- setdiff(observed_species, rownames(trait_data))
  incomplete_traits <- intersect(observed_species, rownames(trait_data)[!stats::complete.cases(trait_data[tt$trait_names])])
  missing_traits <- sort(unique(c(absent_traits, incomplete_traits)))
  if (length(missing_traits) && trait_na == "error") {
    .fd_abort("Species without complete traits: ", paste(missing_traits, collapse = ", "), ".")
  }
  dropped_species <- character()
  if (length(missing_traits)) {
    dropped_species <- missing_traits
    comm <- comm[, setdiff(colnames(comm), missing_traits), drop = FALSE]
  }
  empty_species <- colnames(comm)[colSums(comm) == 0]
  if (length(empty_species)) comm <- comm[, setdiff(colnames(comm), empty_species), drop = FALSE]
  if (!ncol(comm)) .fd_abort("No species with complete traits and positive abundance remain.")
  trait_data <- trait_data[intersect(rownames(trait_data), colnames(comm)), tt$trait_names, drop = FALSE]
  trait_data <- trait_data[colnames(comm), , drop = FALSE]
  trait_observation_info <- .fd_trait_observations(
    trait_observations, tt$trait_names, colnames(comm), itv_na
  )
  structure(list(
    comm = comm,
    meta = meta,
    has_time = isTRUE(has_time),
    traits = trait_data,
    trait_names = tt$trait_names,
    species_col = tt$species_col,
    trait_observations = trait_observation_info$data,
    has_itv = !is.null(trait_observation_info$data),
    itv_na = trait_observation_info$itv_na,
    coord_names = coords,
    latlong = isTRUE(latlong),
    count_data = isTRUE(count_data),
    options = options,
    diagnostics = list(
      missing_traits = missing_traits,
      dropped_species = dropped_species,
      empty_species = empty_species,
      n_trait_observations = if (is.null(trait_observation_info$data)) 0L else nrow(trait_observation_info$data),
      trait_observations_scoped = trait_observation_info$scoped,
      n_samples = nrow(comm),
      n_species = ncol(comm)
    )
  ), class = "fundecomp_input")
}

.fd_zero_itv <- function(species, trait_names) {
  out <- lapply(species, function(.) {
    z <- matrix(0, nrow = 1L, ncol = length(trait_names),
                dimnames = list(NULL, trait_names))
    z
  })
  stats::setNames(out, species)
}

.fd_state_trait_distribution <- function(x, series, time, present_species) {
  centroids <- x$traits
  residuals <- .fd_zero_itv(rownames(centroids), x$trait_names)
  observations <- x$trait_observations
  if (is.null(observations)) return(list(traits = centroids, itv = residuals,
                                         diagnostics = list(fallback_species = character())))
  scoped <- isTRUE(x$diagnostics$trait_observations_scoped)
  if (scoped) {
    observations <- observations[
      observations$series == as.character(series) &
        .fd_time_key(observations$time) == .fd_time_key(time), , drop = FALSE
    ]
  }
  fallback <- character()
  for (sp in present_species) {
    z <- observations[observations$species == sp, x$trait_names, drop = FALSE]
    if (!nrow(z)) {
      if (identical(x$itv_na, "error")) {
        .fd_abort("No trait observations for species `", sp, "` at series `",
                    series, "` and requested time. Use `itv_na = \"species_mean\"` ",
                    "to request an explicit zero-ITV fallback.")
      }
      fallback <- c(fallback, sp)
      next
    }
    centre <- colMeans(z)
    centroids[sp, x$trait_names] <- centre
    residuals[[sp]] <- sweep(as.matrix(z), 2L, centre, "-")
    colnames(residuals[[sp]]) <- x$trait_names
  }
  list(traits = centroids, itv = residuals,
       diagnostics = list(fallback_species = fallback))
}

.fd_state <- function(x, series, time) {
  stopifnot(inherits(x, "fundecomp_input"))
  idx <- x$meta$series == as.character(series) & .fd_time_key(x$meta$time) == .fd_time_key(time)
  if (!any(idx)) .fd_abort("No samples found for series `", series, "` and requested time.")
  comm <- x$comm[idx, , drop = FALSE]
  meta <- x$meta[idx, , drop = FALSE]
  trait_distribution <- .fd_state_trait_distribution(
    x, series, time, colnames(comm)[colSums(comm) > 0]
  )
  list(comm = comm, meta = meta, traits = trait_distribution$traits,
       itv = trait_distribution$itv, itv_diagnostics = trait_distribution$diagnostics,
       trait_names = x$trait_names,
       coord_names = x$coord_names, latlong = x$latlong, count_data = x$count_data,
       series = as.character(series), time = time,
       N = sum(comm), n_samples = nrow(comm), species = colnames(comm))
}

.fd_states <- function(x) {
  keys <- unique(paste(x$meta$series, .fd_time_key(x$meta$time), sep = "@@"))
  first <- match(keys, paste(x$meta$series, .fd_time_key(x$meta$time), sep = "@@"))
  # Preserve the first-appearance order of series, but sort dates within each
  # series by their actual temporal value rather than input row order.
  keys <- unlist(lapply(unique(x$meta$series), function(s) {
    ii <- which(x$meta$series[first] == s)
    keys[ii[.fd_order_time_values(x$meta$time[first[ii]])]]
  }), use.names = FALSE)
  out <- lapply(keys, function(k) {
    z <- strsplit(k, "@@", fixed = TRUE)[[1L]]
    idx <- x$meta$series == z[1L] & .fd_time_key(x$meta$time) == z[2L]
    .fd_state(x, z[1L], x$meta$time[which(idx)[1L]])
  })
  names(out) <- keys
  out
}
