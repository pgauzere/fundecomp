.fd_abort <- function(..., call. = FALSE) {
  stop(paste0(...), call. = call.)
}

# Execute stochastic code without leaking a seed or advancing the caller's
# global RNG stream. This also makes seed = NULL local to the call.
.fd_with_seed <- function(seed, fun) {
  if (!is.function(fun)) .fd_abort("`fun` must be a function.")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.finite(seed)) .fd_abort("`seed` must be one finite value.")
    set.seed(as.integer(seed))
  }
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  force(fun())
}

.fd_effective_seed <- function(seed = NULL) {
  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.finite(seed)) .fd_abort("`seed` must be one finite value.")
    return(as.integer(seed))
  }
  .fd_with_seed(NULL, function() sample.int(.Machine$integer.max, 1L))
}

.fd_check_choice <- function(x, choices, arg) {
  if (length(x) != 1L || is.na(x) || !x %in% choices) {
    .fd_abort("`", arg, "` must be one of: ", paste(choices, collapse = ", "), ".")
  }
  x
}

.fd_as_matrix <- function(x, name = "comm") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    .fd_abort("`", name, "` must be a numeric matrix.")
  }
  if (is.null(rownames(x))) rownames(x) <- paste0("sample", seq_len(nrow(x)))
  if (is.null(colnames(x))) .fd_abort("`", name, "` must have species in its columns.")
  if (anyDuplicated(rownames(x))) .fd_abort("Sample identifiers must be unique.")
  if (anyDuplicated(colnames(x))) .fd_abort("Species identifiers must be unique.")
  x
}

.fd_check_abundance <- function(comm, count_data) {
  if (any(!is.finite(comm)) || any(comm < 0)) {
    .fd_abort("Abundances must be finite and non-negative.")
  }
  if (isTRUE(count_data) && any(comm != floor(comm))) {
    .fd_abort("Count data must contain integer abundances; use `count_data = FALSE` for continuous abundance.")
  }
  invisible(comm)
}

.fd_trait_table <- function(traits, species_col = NULL) {
  if (is.null(traits)) .fd_abort("`traits` is required.")
  traits <- as.data.frame(traits, stringsAsFactors = FALSE)
  if (nrow(traits) == 0L) .fd_abort("`traits` must contain at least one species.")
  if (is.null(species_col)) {
    candidates <- names(traits)[!vapply(traits, is.numeric, logical(1))]
    species_col <- if (length(candidates)) candidates[1L] else names(traits)[1L]
  }
  if (length(species_col) != 1L || !species_col %in% names(traits)) {
    .fd_abort("`species_col` must identify one column in `traits`.")
  }
  sp <- as.character(traits[[species_col]])
  if (anyNA(sp) || any(!nzchar(sp)) || anyDuplicated(sp)) {
    .fd_abort("The species identifier in `traits` must be non-missing and unique.")
  }
  trait_names <- setdiff(names(traits), species_col)
  if (!length(trait_names) || !all(vapply(traits[trait_names], is.numeric, logical(1)))) {
    .fd_abort("`traits` must contain at least one numeric trait column.")
  }
  rownames(traits) <- sp
  traits[[species_col]] <- sp
  list(data = traits, species_col = species_col, trait_names = trait_names)
}

.fd_trait_observations <- function(trait_observations, trait_names, species,
                                     itv_na = c("error", "species_mean")) {
  itv_na <- match.arg(itv_na)
  if (is.null(trait_observations)) {
    return(list(data = NULL, scoped = FALSE, itv_na = itv_na))
  }
  out <- as.data.frame(trait_observations, stringsAsFactors = FALSE)
  required <- c("species", trait_names)
  if (!all(required %in% names(out))) {
    .fd_abort("`trait_observations` must contain: ",
                paste(required, collapse = ", "), ".")
  }
  if (!nrow(out)) .fd_abort("`trait_observations` must contain at least one row.")
  out$species <- as.character(out$species)
  if (anyNA(out$species) || any(!nzchar(out$species))) {
    .fd_abort("Trait-observation species identifiers must be non-missing and non-empty.")
  }
  if (!all(out$species %in% species)) {
    .fd_abort("`trait_observations` contains species absent from the validated community.")
  }
  if (!all(vapply(out[trait_names], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(out[trait_names])))) {
    .fd_abort("Trait-observation columns must be finite numeric values.")
  }
  has_series <- "series" %in% names(out)
  has_time <- "time" %in% names(out)
  if (xor(has_series, has_time)) {
    .fd_abort("Trait observations must contain both `series` and `time`, or neither.")
  }
  if (has_series) {
    out$series <- as.character(out$series)
    if (anyNA(out$series) || any(!nzchar(out$series)) || anyNA(out$time)) {
      .fd_abort("Trait-observation series and time must be non-missing.")
    }
  }
  list(data = out, scoped = has_series, itv_na = itv_na)
}

.fd_time_key <- function(x) {
  # Integer and double representations of the same numeric time must match
  # (for example, metadata created from 1:2 versus a requested time of 1).
  # Keep explicit keys for date-time classes so their printed representation,
  # rather than their storage mode, defines the temporal identity.
  if (inherits(x, "Date")) {
    return(paste0("Date:", as.character(x)))
  }
  if (inherits(x, "POSIXt")) {
    return(paste0("POSIXt:", format(x, tz = "UTC", usetz = TRUE)))
  }
  if (is.numeric(x)) {
    return(paste0("numeric:", format(as.numeric(x), scientific = FALSE,
                                      trim = TRUE, digits = 17)))
  }
  paste0(typeof(x), ":", as.character(x))
}

.fd_order_time_values <- function(x) {
  if (inherits(x, "Date") || inherits(x, "POSIXt") || is.numeric(x)) {
    return(order(x, na.last = NA))
  }
  if (is.ordered(x)) return(order(x, na.last = NA))
  if (is.factor(x) || is.character(x)) return(order(as.character(x), na.last = NA))
  .fd_abort("Unsupported time class; use numeric, Date, POSIXt, factor, or character values.")
}

.fd_state_key <- function(series, time) paste0(as.character(series), "@@", .fd_time_key(time))

.fd_order_times <- function(meta) {
  keys <- unique(.fd_time_key(meta$time))
  keys[order(match(keys, .fd_time_key(meta$time)))]
}

.fd_make_meta <- function(n, series, time, sample, plot_attr = NULL) {
  if (is.null(series)) series <- rep("series1", n)
  if (is.null(time)) time <- rep(1, n)
  if (is.null(sample)) sample <- paste0("sample", seq_len(n))
  if (length(series) != n || length(time) != n || length(sample) != n) {
    .fd_abort("`series`, `time`, and `sample` must have one value per community row.")
  }
  meta <- data.frame(series = as.character(series), time = time, sample = as.character(sample),
                     stringsAsFactors = FALSE)
  if (anyNA(meta$series) || any(!nzchar(meta$series)) || anyNA(meta$time) ||
      anyNA(meta$sample) || any(!nzchar(meta$sample))) {
    .fd_abort("Series, time, and sample identifiers must be non-missing.")
  }
  if (anyDuplicated(paste(meta$series, .fd_time_key(meta$time), meta$sample, sep = "@@"))) {
    .fd_abort("Sample identifiers must be unique within each series and time.")
  }
  if (!is.null(plot_attr)) {
    plot_attr <- as.data.frame(plot_attr, stringsAsFactors = FALSE)
    if (nrow(plot_attr) != n) .fd_abort("`plot_attr` must have one row per community row.")
    keep <- setdiff(names(plot_attr), names(meta))
    meta <- cbind(meta, plot_attr[keep])
  }
  meta
}

.fd_coordinate_names <- function(meta, coord_names) {
  if (is.null(coord_names)) {
    candidates <- intersect(c("x", "y", "X", "Y", "longitude", "latitude", "lon", "lat"), names(meta))
    selected <- if (all(c("longitude", "latitude") %in% candidates)) c("longitude", "latitude")
      else if (all(c("lon", "lat") %in% candidates)) c("lon", "lat")
      else if (all(c("x", "y") %in% candidates)) c("x", "y") else character()
    if (length(selected)) {
      if (!all(vapply(meta[selected], is.numeric, logical(1))) ||
          any(!is.finite(as.matrix(meta[selected])))) {
        .fd_abort("Automatically detected coordinates must be finite numeric columns.")
      }
    }
    return(selected)
  }
  if (length(coord_names) != 2L || !all(coord_names %in% names(meta))) {
    .fd_abort("`coord_names` must name two columns in `plot_attr`.")
  }
  if (!all(vapply(meta[coord_names], is.numeric, logical(1)))) {
    .fd_abort("Coordinate columns must be numeric.")
  }
  if (any(!is.finite(as.matrix(meta[coord_names])))) .fd_abort("Coordinates must be finite.")
  coord_names
}


#' Validate a functional metric decomposition
#'
#' Checks reconstruction, endpoint identity, effort-grid consistency, integer
#' counts, and transition margins that can be verified directly from stored
#' hybrid states. Checks that require identifying a unique TA pairing under
#' ties are reported as not checked; the tie ensemble is retained separately.
#'
#' @param out An `fundecomp_result` returned by [decompose_fmetric_change()].
#' @param tolerance Absolute numerical tolerance for floating-point checks.
#' @return A list with a check table, an overall `passed` flag, maximum error,
#'   and details for checks that could not be evaluated.
#' @export
validate_decomposition <- function(out, tolerance = 1e-8) {
  if (!inherits(out, "fundecomp_result")) .fd_abort("`out` must be an fundecomp_result object.")
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    .fd_abort("`tolerance` must be one finite non-negative number.")
  }
  rows <- list()
  add <- function(check, passed, error = NA_real_, detail = "", status = NULL) {
    rows[[length(rows) + 1L]] <<- data.frame(
      check = check, passed = as.logical(passed), max_error = as.numeric(error),
      detail = detail,
      status = if (is.null(status)) if (is.na(passed)) "unchecked" else
        if (isTRUE(passed)) "passed" else "failed" else status,
      stringsAsFactors = FALSE
    )
  }
  scaled_margin <- function(x, total, count_data) {
    if (sum(x) <= 0 || total <= 0) return(rep(0, length(x)))
    expected <- as.numeric(x) / sum(x) * total
    if (!count_data) return(expected)
    whole <- floor(expected)
    remainder <- as.integer(total - sum(whole))
    if (remainder > 0L) {
      # Largest-remainder apportionment, using the same stable input-order
      # convention as state construction. Comparison below is unlabelled.
      take <- order(-(expected - whole), seq_along(whole))[seq_len(remainder)]
      whole[take] <- whole[take] + 1L
    }
    as.integer(whole)
  }
  margin_error <- function(observed, expected, unlabelled = FALSE) {
    if (length(observed) != length(expected)) return(Inf)
    if (unlabelled) { observed <- sort(as.numeric(observed)); expected <- sort(as.numeric(expected)) }
    if (!length(observed)) return(0)
    max(abs(as.numeric(observed) - as.numeric(expected)))
  }
  # Compare sample profiles within global abundance classes. Species names
  # only align matrices; they do not define equivalence inside tied classes.
  allocation_rank_error <- function(a, b) {
    if (nrow(a$comm) != nrow(b$comm)) return(Inf)
    if (!is.null(rownames(a$comm)) && !is.null(rownames(b$comm))) {
      if (!setequal(rownames(a$comm), rownames(b$comm))) return(Inf)
      b$comm <- b$comm[rownames(a$comm), , drop = FALSE]
    } else if (is.null(rownames(a$comm)) != is.null(rownames(b$comm))) {
      return(Inf)
    }
    ca <- colSums(a$comm); cb <- colSums(b$comm)
    if (length(ca) != length(cb)) return(Inf)
    if (max(abs(sort(ca) - sort(cb))) > tolerance) return(Inf)
    # Construct abundance classes with the validator tolerance. Counterfactual
    # rescaling can represent an intended tie (e.g. 4.2) as values differing
    # by a few floating-point ulps. Exact `unique()`/`==` would split that tie
    # and falsely fail the TA allocation check.
    abundance_classes <- function(x) {
      ord <- order(x)
      sorted <- x[ord]
      if (!length(sorted)) return(list())
      groups <- list()
      current <- ord[1L]
      anchor <- sorted[1L]
      if (length(sorted) > 1L) for (j in seq.int(2L, length(sorted))) {
        if (abs(sorted[j] - anchor) <= tolerance) {
          current <- c(current, ord[j])
        } else {
          groups[[length(groups) + 1L]] <- current
          current <- ord[j]
          anchor <- sorted[j]
        }
      }
      groups[[length(groups) + 1L]] <- current
      groups
    }
    classes_a <- abundance_classes(ca)
    classes_b <- abundance_classes(cb)
    if (length(classes_a) != length(classes_b)) return(Inf)
    profile_match_error <- function(x, y) {
      if (!identical(dim(x), dim(y))) return(Inf)
      if (!ncol(x)) return(0)
      # Compare profiles without relying on a canonical sort (whose ties can
      # inherit column order). A maximum bipartite match identifies whether
      # the two classes contain the same multiset of sample profiles.
      distance <- t(vapply(seq_len(ncol(x)), function(i) {
        vapply(seq_len(ncol(y)), function(j) {
          den <- pmax(1, abs(x[, i]), abs(y[, j]))
          max(abs(x[, i] - y[, j]) / den)
        }, numeric(1))
      }, numeric(ncol(y))))
      adjacency <- lapply(seq_len(ncol(x)), function(i)
        which(distance[i, ] <= tolerance))
      matched_x <- integer(ncol(y))
      augment <- function(i, seen) {
        for (j in adjacency[[i]]) {
          if (seen[j]) next
          seen[j] <- TRUE
          if (matched_x[j] == 0L || augment(matched_x[j], seen)) {
            matched_x[j] <<- i
            return(TRUE)
          }
        }
        FALSE
      }
      for (i in seq_len(ncol(x))) if (!augment(i, rep(FALSE, ncol(y)))) return(Inf)
      paired_y <- which(matched_x > 0L)
      max(distance[cbind(matched_x[paired_y], paired_y)])
    }
    errs <- numeric()
    for (i in seq_along(classes_a)) {
      ia <- classes_a[[i]]; ib <- classes_b[[i]]
      if (length(ia) != length(ib)) return(Inf)
      # Match by abundance value rather than class ordinal, which may depend on
      # the input column ordering.
      va <- mean(ca[ia]); vb <- mean(cb[ib])
      if (abs(va - vb) > tolerance) return(Inf)
      xa <- a$comm[, ia, drop = FALSE]
      xb <- b$comm[, ib, drop = FALSE]
      errs <- c(errs, profile_match_error(xa, xb))
    }
    if (length(errs)) max(errs) else 0
  }
  trait_margin_error <- function(a, b) {
    ia <- colSums(a$comm) > 0; ib <- colSums(b$comm) > 0
    xa <- as.matrix(a$traits[ia, , drop = FALSE])
    xb <- as.matrix(b$traits[ib, , drop = FALSE])
    if (!identical(dim(xa), dim(xb))) return(NA_real_)
    ord_rows <- function(z) {
      if (nrow(z) < 2L) return(z)
      z[do.call(order, c(unname(as.list(as.data.frame(z))), list(method = "radix"))), , drop = FALSE]
    }
    xa <- ord_rows(xa); xb <- ord_rows(xb)
    max(abs(xa - xb), na.rm = TRUE)
  }
  itv_margin_key <- function(state) {
    spp <- names(which(colSums(state$comm) > 0))
    keys <- vapply(spp, function(sp) {
      z <- as.matrix(state$itv[[sp]])
      if (nrow(z) > 1L) z <- z[do.call(order, c(unname(as.list(as.data.frame(z))),
                                                  list(method = "radix"))), , drop = FALSE]
      paste(nrow(z), ncol(z), paste(format(as.numeric(z), digits = 16,
        scientific = TRUE), collapse = ","), sep = "|")
    }, character(1))
    sort(keys)
  }
  recon_values <- if ("reconstruction_error" %in% names(out$curves))
    abs(out$curves$reconstruction_error) else numeric()
  recon_values <- recon_values[is.finite(recon_values)]
  recon <- if (length(recon_values)) max(recon_values) else NA_real_
  add("reconstruction_identity", is.finite(recon) && recon <= tolerance,
      recon, if (isTRUE(out$full_decomposition))
        "Contributions reconstruct the observed endpoint contrast." else
        "Selected contributions reconstruct the conditional partial contrast.")

  if (!is.null(out$states) && length(out$states)) {
    endpoint_error <- function(a, b) {
      if (!identical(dim(a), dim(b)) || !identical(dimnames(a), dimnames(b))) return(Inf)
      max(abs(a - b), na.rm = TRUE)
    }
    e0 <- endpoint_error(out$states[[1L]]$comm, out$endpoint$reference$comm)
    full_decomposition <- isTRUE(out$full_decomposition)
    if (full_decomposition) {
      e1 <- endpoint_error(out$states[[length(out$states)]]$comm, out$endpoint$focal$comm)
    } else {
      e1 <- NA_real_
    }
    add("reference_endpoint", is.finite(e0) && e0 <= tolerance, e0,
        "First stored path state equals the observed reference community.")
    add(if (full_decomposition) "focal_endpoint" else "partial_counterfactual_endpoint",
        if (full_decomposition) is.finite(e1) && e1 <= tolerance else NA,
        e1, if (full_decomposition) "Last stored path state equals the observed focal community." else
          "Last stored path state is the selected-components counterfactual; it is not asserted to equal the observed focal community.")
    trait_error <- function(a, b) {
      a <- as.matrix(a); b <- as.matrix(b)
      if (!identical(dim(a), dim(b)) || !identical(dimnames(a), dimnames(b))) return(Inf)
      max(abs(a - b), na.rm = TRUE)
    }
    t0 <- trait_error(out$states[[1L]]$traits, out$endpoint$reference$traits)
    t1 <- if (full_decomposition) trait_error(out$states[[length(out$states)]]$traits,
      out$endpoint$focal$traits) else NA_real_
    add("reference_endpoint_STD", is.finite(t0) && t0 <= tolerance, t0,
        "First path state retains the observed reference trait rows.")
    add(if (full_decomposition) "focal_endpoint_STD" else "partial_counterfactual_STD",
        if (full_decomposition) is.finite(t1) && t1 <= tolerance else NA,
        t1, if (full_decomposition) "Last path state retains the observed focal trait rows." else
          "The final counterfactual trait rows need not equal the observed focal trait rows when STD or ITV is held at reference.")
    states <- out$states
    if (length(states)) {
      integer_ok <- all(vapply(states, function(s) !isTRUE(s$count_data) ||
        all(is.finite(s$comm) & s$comm == floor(s$comm)), logical(1)))
      add("integer_counts", integer_ok, if (integer_ok) 0 else Inf,
          "All stored count states contain integer abundances.")
    }
  } else {
    add("stored_path_states", NA, NA_real_, "Path states are unavailable for this Shapley result.")
  }
  effort_ok <- length(unique(out$curves$effort)) == length(out$effort) &&
    setequal(unique(out$curves$effort), out$effort)
  add("common_effort_grid", effort_ok, if (effort_ok) 0 else Inf,
      "All components use the declared common effort grid.")

  if (!is.null(out$states) && length(out$states) > 1L) {
    comp <- out$path
    for (j in seq_along(comp)) {
      before <- out$states[[j]]
      after <- out$states[[j + 1L]]
      component <- comp[j]
      # Total abundance and global SAD are checked for every transition.
      nerr <- if (component == "N") abs(sum(after$comm) - sum(out$endpoint$focal$comm)) else
        abs(sum(after$comm) - sum(before$comm))
      add(paste0(component, "_N"), nerr <= tolerance, nerr,
          if (component == "N") "N changes toward the observed focal total." else "N is held fixed.")
      c0 <- colSums(before$comm); c1 <- colSums(after$comm)
      saderr <- if (component == "SAD") NA_real_ else if (component == "N") {
        expected <- scaled_margin(c0, sum(after$comm),
                                  isTRUE(before$count_data) && isTRUE(after$count_data))
        margin_error(c1, expected, unlabelled = TRUE)
      } else margin_error(c1, c0)
      add(paste0(component, "_SAD"), if (is.na(saderr)) NA else saderr <= tolerance,
          saderr, if (component == "SAD") "SAD changes at this transition." else if (component == "N")
            "Unlabelled SAD is compared with the deterministically scaled reference margin." else
            "Global species-abundance margins are held fixed.")
      if (component != "aggregation") {
        r0 <- rowSums(before$comm); r1 <- rowSums(after$comm)
        if (component == "N") {
          expected_rows <- scaled_margin(r0, sum(after$comm),
                                         isTRUE(before$count_data) && isTRUE(after$count_data))
          aerr <- margin_error(r1, expected_rows)
        } else aerr <- if (length(r0) == length(r1)) max(abs(r0 - r1), na.rm = TRUE) else Inf
        add(paste0(component, "_sample_allocation"), aerr <= tolerance, aerr,
            if (component == "N") "Sample totals match deterministic rescaling of reference row margins." else
              "Sample-total profile is held fixed.")
      } else {
        add(paste0(component, "_sample_allocation"), NA, NA_real_,
            "Aggregation may change sample totals and sample-species allocation.")
      }
      if (component != "traits") {
        terr <- trait_margin_error(before, after)
        add(paste0(component, "_STD"), if (is.na(terr)) NA else terr <= tolerance,
            terr, if (is.na(terr))
              "STD row counts differ; systematic transport across richness is approximate." else
              "Unlabelled complete trait rows are held fixed, checked without species labels.")
      }
      if ("ITV" %in% out$path) {
        if (component != "ITV") {
          ka <- itv_margin_key(before); kb <- itv_margin_key(after)
          if (component == "SAD" && sum(colSums(before$comm) > 0) != sum(colSums(after$comm) > 0)) {
            add(paste0(component, "_ITV"), NA, NA_real_,
                "Positive richness changes, so equality of transported ITV margins is not defined.",
                status = "approximate_systematic_transport")
          } else {
            same_itv <- identical(ka, kb)
            add(paste0(component, "_ITV"), same_itv, if (same_itv) 0 else Inf,
                "Unlabelled within-species residual-vector distributions are compared without species labels.")
          }
        }
      }
      previous_sources <- before$construction$sources
      if (is.null(previous_sources)) previous_sources <- stats::setNames(
        rep("reference", 6L), c("N", "aggregation", "SAD", "STD", "ITV", "TA"))
      current_sources <- after$construction$sources
      if (!is.null(current_sources) && "TA" %in% names(current_sources)) {
        old_ta <- unname(previous_sources["TA"])
        new_ta <- unname(current_sources["TA"])
        expected_ta <- if (component == "TA") "focal" else old_ta
        ta_rule_ok <- identical(new_ta, expected_ta)
        add(paste0(component, "_TA_rule"), ta_rule_ok,
            if (ta_rule_ok) 0 else Inf,
            "Checks the recorded source of the TA replacement rule; tied estimates are averaged over pairings.")
      }
      if (component %in% c("traits", "ITV")) {
        commerr <- if (identical(dim(before$comm), dim(after$comm)) &&
                       identical(dimnames(before$comm), dimnames(after$comm)))
          max(abs(before$comm - after$comm), na.rm = TRUE) else Inf
        add(paste0(component, "_community_fixed"), commerr <= tolerance, commerr,
            "Species-labelled community is fixed for STD and ITV transitions.")
      }
      if (component == "TA") {
        commerr <- allocation_rank_error(before, after)
        add("TA_community_fixed", commerr <= tolerance, commerr,
            "N, unlabelled SAD, sample totals, and sample profiles within abundance-rank classes are preserved modulo species relabelling.")
      }
    }
  }
  check <- if (length(rows)) do.call(rbind, rows) else data.frame()
  tested <- !is.na(check$passed)
  passed <- all(check$passed[tested]) && any(tested) &&
    all(is.na(check$passed) | check$passed)
  errors <- check$max_error[is.finite(check$max_error)]
  structure(list(check = check, passed = passed,
                 max_error = if (length(errors)) max(errors) else NA_real_,
                 details = list(n_checked = sum(tested), n_not_checked = sum(!tested),
                                tolerance = tolerance)), class = "fundecomp_validation")
}

.fd_transition_constraint_table <- function(out) {
  path <- out$path
  properties <- c("N", "sample-level allocation", "SAD", "STD", "ITV", "TA")
  if (is.null(path)) return(data.frame())
  rows <- lapply(path, function(component) {
    changed <- switch(component, N = "N", aggregation = "sample-level allocation",
                      SAD = "SAD", traits = "STD", ITV = "ITV", TA = "TA")
    action <- ifelse(properties == changed, "changed_toward_focal", "held_fixed")
    if (component == "aggregation") {
      action[properties == "sample-level allocation"] <- "changed_toward_focal"
      action[properties == "N"] <- "conditioned_upon"
      action[properties == "SAD"] <- "conditioned_upon"
    }
    if (component == "SAD") action[properties == "sample-level allocation"] <- "held_fixed"
    if (!"ITV" %in% path) action[properties == "ITV"] <- "not_applicable"
    data.frame(transition = component,
      property = properties, action = action,
      preservation = ifelse(action == "conditioned_upon", "conditioned",
        ifelse(action == "not_applicable", "not_applicable",
        ifelse(component == "aggregation" & properties == "sample-level allocation",
               "may_change_by_definition",
        ifelse(component == "SAD" & properties %in% c("STD", "ITV"),
               "approximate_systematic_transport", "direct_check_if_available")))),
      checked = FALSE, passed = NA, status = "unchecked", max_difference = NA_real_,
      note = ifelse(properties == "sample-level allocation" & component == "aggregation",
        "Includes changes in sample totals and species allocation.",
        ifelse(properties == "TA",
          "Validation checks the recorded TA replacement-rule source; pairing identity is not inferred from one tied representative state.",
        ifelse(properties == "ITV" & !("ITV" %in% path),
          "ITV is absent because individual trait observations were not supplied.",
        ifelse(component == "SAD" & properties %in% c("STD", "ITV"),
          "Trait margins are systematically transported if positive richness changes; this is not an exact identity check.",
          "See validate_decomposition(out) for programmatic checks.")))),
      stringsAsFactors = FALSE)
  })
  table <- do.call(rbind, rows)
  validation <- validate_decomposition(out)
  # Attach checks by their explicit transition/property names where a direct
  # check exists. Unchecked properties remain visibly unchecked.
  for (i in seq_len(nrow(table))) {
    transition_key <- table$transition[i]
    property_key <- switch(table$property[i],
      `sample-level allocation` = "sample_allocation",
      TA = "TA_rule", table$property[i])
    key <- paste0(transition_key, "_", property_key)
    if (transition_key == "TA" && table$property[i] == "sample-level allocation") {
      key <- "TA_community_fixed"
    }
    hit <- match(key, validation$check$check)
    if (!is.na(hit)) {
      table$checked[i] <- !is.na(validation$check$passed[hit])
      table$passed[i] <- validation$check$passed[hit]
      table$status[i] <- validation$check$status[hit]
      table$max_difference[i] <- validation$check$max_error[hit]
    }
  }
  table
}
