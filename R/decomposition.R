.fd_components <- c("N", "aggregation", "SAD", "traits", "TA")
.fd_components_with_itv <- c("N", "aggregation", "SAD", "traits", "ITV", "TA")
.fd_component_labels <- c(N = "N", aggregation = "Aggregation", SAD = "SAD",
                            traits = "STD", ITV = "ITV", TA = "Trait-Abundance")

.fd_components_for <- function(x) {
  if (isTRUE(x$has_itv)) .fd_components_with_itv else .fd_components
}

.fd_union_state <- function(state, species, traits) {
  comm <- matrix(0, nrow = nrow(state$comm), ncol = length(species),
                 dimnames = list(rownames(state$comm), species))
  common <- intersect(colnames(state$comm), species)
  comm[, common] <- state$comm[, common, drop = FALSE]
  itv <- .fd_zero_itv(species, names(traits))
  common_itv <- intersect(names(state$itv), species)
  if (length(common_itv)) itv[common_itv] <- state$itv[common_itv]
  list(comm = comm, meta = state$meta, traits = traits[species, , drop = FALSE], itv = itv,
       trait_names = state$trait_names, coord_names = state$coord_names,
       latlong = state$latlong, count_data = isTRUE(state$count_data),
       series = state$series, time = state$time,
       N = sum(comm), n_samples = nrow(comm), species = species)
}

.fd_state_traits <- function(x, species) {
  out <- x$traits[match(species, rownames(x$traits)), , drop = FALSE]
  out <- out[species, , drop = FALSE]
  if (anyNA(out)) .fd_abort("Trait data are incomplete for the comparison species.")
  attr(out, "species_col") <- x$species_col
  out
}

.fd_rescale_rows <- function(comm, N) {
  if (sum(comm) <= 0 || N <= 0) return(comm * 0)
  comm * (N / sum(comm))
}

.fd_ipf <- function(pattern, row_targets, col_targets, maxit = 5000L, tol = 1e-12) {
  z <- as.matrix(pattern)
  z[!is.finite(z) | z < 0] <- 0
  z <- z + 1e-12
  for (i in seq_len(maxit)) {
    row_factor <- ifelse(rowSums(z) > 0, row_targets / rowSums(z), 0)
    row_factor[!is.finite(row_factor)] <- 0
    z <- sweep(z, 1L, as.numeric(row_factor), "*")
    col_factor <- ifelse(colSums(z) > 0, col_targets / colSums(z), 0)
    col_factor[!is.finite(col_factor)] <- 0
    z <- sweep(z, 2L, as.numeric(col_factor), "*")
    if (max(abs(rowSums(z) - row_targets)) < tol && max(abs(colSums(z) - col_targets)) < tol) break
  }
  z
}

.fd_integer_margin <- function(weights, total) {
  if (length(total) != 1L || !is.finite(total) || total < 0 || total != floor(total)) {
    .fd_abort("Counterfactual total abundance must be a non-negative integer.")
  }
  weights <- as.numeric(weights)
  weights[!is.finite(weights) | weights < 0] <- 0
  if (!length(weights)) return(integer())
  if (total == 0L) return(rep.int(0L, length(weights)))
  if (sum(weights) <= 0) weights <- rep(1, length(weights))
  expected <- weights / sum(weights) * total
  out <- floor(expected)
  remainder <- as.integer(total - sum(out))
  if (remainder > 0L) {
    take <- order(-(expected - out), seq_along(out))[seq_len(remainder)]
    out[take] <- out[take] + 1L
  }
  as.integer(out)
}

.fd_integerize_ipf <- function(expected, row_targets, col_targets) {
  row_targets <- as.integer(row_targets)
  col_targets <- as.integer(col_targets)
  if (sum(row_targets) != sum(col_targets)) {
    .fd_abort("Integer row and species margins do not have the same total.")
  }
  out <- floor(expected)
  storage.mode(out) <- "integer"
  row_deficit <- row_targets - rowSums(out)
  col_deficit <- col_targets - colSums(out)
  if (any(row_deficit < 0L) || any(col_deficit < 0L)) {
    .fd_abort(
      "IPF rounding produced an invalid negative margin deficit ",
      "(minimum row deficit ", min(row_deficit),
      ", minimum column deficit ", min(col_deficit),
      "; maximum pre-rounding margin error ",
      signif(max(abs(c(rowSums(expected) - row_targets,
                       colSums(expected) - col_targets))), 6), ")."
    )
  }
  fractional <- expected - floor(expected)
  additions <- matrix(0, nrow(expected), ncol(expected))
  while (sum(row_deficit) > 0L) {
    eligible <- outer(row_deficit > 0L, col_deficit > 0L, "&")
    if (!any(eligible)) .fd_abort("Could not integerize the IPF matrix with exact margins.")
    score <- fractional - additions
    score[!eligible] <- -Inf
    cell <- arrayInd(which.max(score), dim(score))
    out[cell[1L], cell[2L]] <- out[cell[1L], cell[2L]] + 1L
    additions[cell[1L], cell[2L]] <- additions[cell[1L], cell[2L]] + 1
    row_deficit[cell[1L]] <- row_deficit[cell[1L]] - 1L
    col_deficit[cell[2L]] <- col_deficit[cell[2L]] - 1L
  }
  dimnames(out) <- dimnames(expected)
  if (!identical(as.integer(rowSums(out)), row_targets) ||
      !identical(as.integer(colSums(out)), col_targets)) {
    .fd_abort("Integerized counterfactual state does not preserve its target margins.")
  }
  out
}

.fd_integer_counterfactual <- function(pattern, row_weights, col_weights, N) {
  row_targets <- .fd_integer_margin(row_weights, N)
  col_targets <- .fd_integer_margin(col_weights, N)
  expected <- .fd_ipf(pattern, row_targets, col_targets)
  out <- .fd_integerize_ipf(expected, row_targets, col_targets)
  attr(out, "construction") <- list(
    method = "integerized_ipf",
    row_targets = stats::setNames(row_targets, rownames(pattern)),
    col_targets = stats::setNames(col_targets, colnames(pattern)),
    max_expected_rounding_error = max(abs(out - expected))
  )
  out
}

.fd_continuous_counterfactual <- function(pattern, row_weights, col_weights, N) {
  row_weights <- as.numeric(row_weights)
  col_weights <- as.numeric(col_weights)
  if (sum(row_weights) <= 0 || sum(col_weights) <= 0 || N < 0) {
    .fd_abort("Continuous counterfactual margins must have positive weights and a non-negative total.")
  }
  row_targets <- row_weights / sum(row_weights) * N
  col_targets <- col_weights / sum(col_weights) * N
  out <- .fd_ipf(pattern, row_targets, col_targets)
  attr(out, "construction") <- list(
    method = "continuous_ipf",
    row_targets = stats::setNames(row_targets, rownames(pattern)),
    col_targets = stats::setNames(col_targets, colnames(pattern)),
    max_margin_error = max(abs(c(rowSums(out) - row_targets,
                                 colSums(out) - col_targets)))
  )
  out
}

.fd_systematic_trait_indices <- function(n, k, weights = NULL) {
  if (n < 1L || k < 1L) return(integer())
  if (is.null(weights)) weights <- rep(1 / n, n)
  weights <- as.numeric(weights)
  weights[!is.finite(weights) | weights < 0] <- 0
  if (sum(weights) <= 0) weights <- rep(1 / n, n) else weights <- weights / sum(weights)
  u <- (seq_len(k) - 0.5) / k
  pmax(1L, pmin(n, findInterval(u, cumsum(weights)) + 1L))
}

.fd_rank_order <- function(abundance) {
  abundance <- abundance[is.finite(abundance) & abundance > 0]
  if (!length(abundance)) return(character())
  names(abundance)[order(-abundance, names(abundance), method = "radix")]
}

# Abundance ranks are only identifiable between distinct abundance values.
# These helpers keep the historical order available for strict SADs, while
# representing tied ranks as the set of all admissible within-class orders.
.fd_tie_components <- function(abundance, key = NULL) {
  nms <- names(abundance)
  if (is.null(nms)) .fd_abort("Abundance vectors used for tie handling must be named.")
  abundance <- stats::setNames(as.numeric(abundance), nms)
  positive <- is.finite(abundance) & abundance > 0
  if (!any(positive)) return(list(base = character(), groups = list(), count = 1))
  base <- .fd_rank_order(abundance)
  split_id <- match(as.character(abundance[base]), unique(as.character(abundance[base])))
  groups <- split(base, split_id)

  # A data-derived key is used only to make Monte Carlo draws stable under
  # relabelling. It is not a scientific ranking: all within-class
  # permutations are still admissible and are averaged or sampled uniformly.
  if (!is.null(key)) {
    key <- as.data.frame(key, stringsAsFactors = FALSE)
    if (is.null(rownames(key)) || !all(base %in% rownames(key))) {
      .fd_abort("Tie-breaking keys must have rows named by the abundance vector.")
    }
    key <- key[base, , drop = FALSE]
    groups <- lapply(groups, function(group) {
      values <- key[group, , drop = FALSE]
      ord <- do.call(order, c(unname(as.list(values)), list(seq_along(group)),
                              list(method = "radix")))
      group[ord]
    })
  }
  count <- 1
  for (group in groups) {
    size <- length(group)
    if (size > 1L) {
      count <- if (count > .Machine$integer.max / factorial(size)) Inf else
        count * factorial(size)
    }
  }
  list(base = unlist(groups, use.names = FALSE), groups = groups, count = count)
}

.fd_permute_group <- function(group) {
  if (length(group) <= 1L) return(list(group))
  out <- list()
  for (i in seq_along(group)) {
    rest <- group[-i]
    for (tail in .fd_permute_group(rest)) {
      out[[length(out) + 1L]] <- c(group[i], tail)
    }
  }
  out
}

.fd_all_tie_orders <- function(groups) {
  orders <- list(character())
  for (group in groups) {
    group_orders <- .fd_permute_group(group)
    orders <- unlist(lapply(orders, function(prefix) {
      lapply(group_orders, function(suffix) c(prefix, suffix))
    }), recursive = FALSE)
  }
  orders
}

.fd_sample_tie_order <- function(groups) {
  unlist(lapply(groups, function(group) {
    if (length(group) <= 1L) group else sample(group, length(group), replace = FALSE)
  }), use.names = FALSE)
}

.fd_tie_plan <- function(source_abundance, target_abundance, source_key = NULL,
                           target_key = NULL, tie_method = c("average", "random", "legacy"),
                           n_tie_perm = 200L, max_exact_tie_perm = 10000L,
                           seed = NULL) {
  tie_method <- match.arg(tie_method)
  if (tie_method == "legacy") {
    source_key <- NULL
    target_key <- NULL
  }
  if (length(n_tie_perm) != 1L || !is.finite(n_tie_perm) || n_tie_perm < 1L ||
      n_tie_perm != floor(n_tie_perm)) {
    .fd_abort("`n_tie_perm` must be one positive integer.")
  }
  if (length(max_exact_tie_perm) != 1L || !is.finite(max_exact_tie_perm) ||
      max_exact_tie_perm < 1L || max_exact_tie_perm != floor(max_exact_tie_perm)) {
    .fd_abort("`max_exact_tie_perm` must be one positive integer.")
  }
  source <- .fd_tie_components(source_abundance, source_key)
  target <- .fd_tie_components(target_abundance, target_key)
  total <- if (is.infinite(source$count) || is.infinite(target$count) ||
               source$count > max_exact_tie_perm / target$count) Inf else
    source$count * target$count
  has_ties <- any(vapply(c(source$groups, target$groups), length, integer(1)) > 1L)

  if (tie_method == "legacy" || !has_ties) {
    orders <- list(list(source = source$base, target = target$base))
    return(list(orders = orders, total = total, exact = !has_ties,
                method = if (has_ties) tie_method else "none",
                n_used = 1L, tied = has_ties,
                source_group_sizes = unlist(lapply(source$groups, length),
                                            use.names = FALSE),
                target_group_sizes = unlist(lapply(target$groups, length),
                                            use.names = FALSE),
                group_sizes = unlist(c(lapply(source$groups, length),
                                       lapply(target$groups, length)),
                                     use.names = FALSE)))
  }

  exact <- identical(tie_method, "average") && is.finite(total) &&
    total <= max_exact_tie_perm
  orders <- .fd_with_seed(seed, function() {
    if (exact) {
      source_orders <- .fd_all_tie_orders(source$groups)
      target_orders <- .fd_all_tie_orders(target$groups)
      out <- vector("list", length(source_orders) * length(target_orders))
      k <- 0L
      for (so in source_orders) for (to in target_orders) {
        k <- k + 1L
        out[[k]] <- list(source = so, target = to)
      }
      out
    } else {
      n <- if (tie_method == "random") 1L else as.integer(n_tie_perm)
      replicate(n, list(source = .fd_sample_tie_order(source$groups),
                        target = .fd_sample_tie_order(target$groups)),
                simplify = FALSE)
    }
  })
  list(orders = orders, total = total, exact = exact, method = tie_method,
       n_used = length(orders), tied = has_ties,
       source_group_sizes = unlist(lapply(source$groups, length),
                                   use.names = FALSE),
       target_group_sizes = unlist(lapply(target$groups, length),
                                   use.names = FALSE),
       group_sizes = unlist(c(lapply(source$groups, length),
                              lapply(target$groups, length)),
                            use.names = FALSE))
}

.fd_ta_axis_order <- function(species, ta_axis_scores, axis_order = NULL) {
  if (!length(species)) return(character())
  if (is.null(ta_axis_scores) || is.null(names(ta_axis_scores)) ||
      !all(species %in% names(ta_axis_scores))) {
    .fd_abort("A validated scalar T-A axis is required to order complete trait rows.")
  }
  scores <- as.numeric(ta_axis_scores[species])
  if (any(!is.finite(scores))) {
    .fd_abort("T-A axis scores must be finite for every transported species.")
  }
  if (!is.null(axis_order)) {
    if (!all(species %in% axis_order)) .fd_abort("The tie-averaged TA-axis order omits transported species.")
    return(axis_order[axis_order %in% species])
  }
  # This fallback is used only for observed endpoints and internal prototypes.
  # Counterfactual states average every within-score ordering in
  # .fd_subset_state_variants().
  species[order(scores, method = "radix")]
}

.fd_axis_tie_components <- function(scores, traits) {
  species <- names(scores)
  if (is.null(species) || !length(species) || anyDuplicated(species)) {
    .fd_abort("Trait-axis tie scores must be named by unique species.")
  }
  scores <- scores[species]
  if (any(!is.finite(scores))) .fd_abort("Trait-axis tie scores must be finite.")
  ord <- order(scores, method = "radix")
  ordered <- species[ord]
  groups <- split(ordered, match(as.character(scores[ordered]), unique(as.character(scores[ordered]))))
  traits <- as.data.frame(traits[species, , drop = FALSE])
  # Data-derived ordering is only a stable enumeration convention. Every
  # within-axis-tie ordering is still included in exact averaging or sampled.
  groups <- lapply(groups, function(group) {
    if (length(group) < 2L) return(group)
    values <- traits[group, , drop = FALSE]
    group[do.call(order, c(unname(as.list(values)), list(method = "radix")))]
  })
  count <- 1
  for (group in groups) {
    if (length(group) > 1L) {
      f <- factorial(length(group))
      count <- if (is.infinite(count) || count > .Machine$integer.max / f) Inf else count * f
    }
  }
  list(base = unlist(groups, use.names = FALSE), groups = groups,
       group_sizes = lengths(groups), count = count,
       tied = any(lengths(groups) > 1L))
}

.fd_axis_tie_orders <- function(plan, exact, n, seed) {
  if (!plan$tied) return(list(plan$base))
  .fd_with_seed(seed, function() {
    if (exact) return(.fd_all_tie_orders(plan$groups))
    replicate(n, .fd_sample_tie_order(plan$groups), simplify = FALSE)
  })
}

.fd_species_tie_key <- function(state, scores, traits) {
  spp <- names(scores)
  z <- as.data.frame(traits[spp, , drop = FALSE], stringsAsFactors = FALSE)
  names(z) <- paste0("trait_", seq_len(ncol(z)))
  z$ta_axis <- as.numeric(scores[spp])
  profiles <- t(vapply(spp, function(sp) sort(as.numeric(state$comm[, sp])), numeric(nrow(state$comm))))
  profiles <- as.data.frame(profiles, stringsAsFactors = FALSE)
  names(profiles) <- paste0("sample_margin_", seq_len(ncol(profiles)))
  cbind(z, profiles)
}

.fd_abundance_slots <- function(source, species, N, count_data) {
  raw <- colSums(source$comm)
  ranked_species <- .fd_rank_order(raw)
  if (!length(ranked_species)) {
    .fd_abort("A counterfactual SAD source contains no individuals.")
  }
  margin <- if (count_data) {
    .fd_integer_margin(raw[ranked_species], N)
  } else {
    raw[ranked_species] / sum(raw[ranked_species]) * N
  }
  out <- stats::setNames(rep(0, length(species)), species)
  out[ranked_species] <- margin
  out
}

.fd_ranked_aggregation_pattern <- function(source, target_abundance) {
  source_abundance <- colSums(source$comm)
  source_species <- .fd_rank_order(source_abundance)
  target_species <- .fd_rank_order(target_abundance)
  if (!length(source_species) || !length(target_species)) {
    .fd_abort("Aggregation patterns require positive source and target abundance.")
  }
  profiles <- sweep(
    source$comm[, source_species, drop = FALSE], 2L,
    source_abundance[source_species], "/"
  )
  indices <- .fd_systematic_trait_indices(
    length(source_species), length(target_species)
  )
  out <- matrix(
    0, nrow = nrow(source$comm), ncol = length(target_abundance),
    dimnames = list(rownames(source$comm), names(target_abundance))
  )
  out[, target_species] <- profiles[, indices, drop = FALSE]
  out
}

.fd_pairing_profile <- function(source, target_richness,
                                  abundance_order = NULL,
                                  ta_axis_scores = NULL,
                                  axis_order = NULL) {
  source_abundance <- colSums(source$comm)
  if (is.null(abundance_order)) abundance_order <- .fd_rank_order(source_abundance)
  if (!length(abundance_order) || target_richness < 1L) return(numeric())
  axis_order <- .fd_ta_axis_order(abundance_order, ta_axis_scores, axis_order)
  # This permutation is the copula-like association between abundance ranks
  # and positions on the declared scalar T-A axis. Complete multivariate rows
  # are transported using this order below; the axis is not used as the metric
  # trait space.
  trait_positions <- match(abundance_order, axis_order)
  profile <- (trait_positions - 0.5) / length(trait_positions)
  indices <- .fd_systematic_trait_indices(length(profile), target_richness)
  profile[indices]
}

.fd_ranked_trait_margin <- function(std_source, pairing_source, target_abundance,
                                      pairing_order = NULL, target_order = NULL,
                                      ta_axis_scores = NULL, axis_order = NULL) {
  traits <- std_source$traits
  if (is.null(target_order)) target_order <- .fd_rank_order(target_abundance)
  target_species <- target_order
  std_species <- .fd_rank_order(colSums(std_source$comm))
  if (!length(target_species) || !length(std_species)) {
    .fd_abort("STD construction requires positive source and target richness.")
  }

  std_axis_order <- .fd_ta_axis_order(std_species, ta_axis_scores, axis_order)
  margin_indices <- .fd_systematic_trait_indices(
    length(std_axis_order), length(target_species)
  )
  selected_margin <- traits[std_axis_order[margin_indices], , drop = FALSE]
  if (is.null(pairing_order)) pairing_order <- .fd_rank_order(colSums(pairing_source$comm))
  pairing_profile <- .fd_pairing_profile(
    pairing_source, length(target_species), pairing_order, ta_axis_scores,
    axis_order
  )
  pairing_permutation <- rank(pairing_profile, ties.method = "first")

  out <- traits
  out[target_species, ] <- selected_margin[pairing_permutation, , drop = FALSE]
  rownames(out) <- rownames(traits)
  attr(out, "species_col") <- attr(traits, "species_col")
  out
}

.fd_ranked_itv_margin <- function(itv_source, target_traits, target_abundance,
                                    target_order = NULL, ta_axis_scores = NULL,
                                    axis_order = NULL) {
  if (is.null(target_order)) target_order <- .fd_rank_order(target_abundance)
  target_species <- target_order
  source_species <- .fd_rank_order(colSums(itv_source$comm))
  out <- .fd_zero_itv(rownames(target_traits), names(target_traits))
  if (!length(target_species) || !length(source_species)) return(out)
  source_order <- .fd_ta_axis_order(source_species, ta_axis_scores, axis_order)
  target_order <- .fd_ta_axis_order(target_species, ta_axis_scores, axis_order)
  indices <- .fd_systematic_trait_indices(length(source_order), length(target_order))
  out[target_order] <- itv_source$itv[source_order[indices]]
  out
}

.fd_resolve_ta_axis <- function(ta_axis, traits, trait_names,
                                  analysed_trait = NULL) {
  trait_names <- as.character(trait_names)
  if (!length(trait_names)) {
    .fd_abort("A T-A axis cannot be defined without trait columns.")
  }
  selected <- analysed_trait
  if (is.null(selected)) {
    selected <- trait_names
  } else if (is.numeric(selected)) {
    if (any(!is.finite(selected)) || any(selected != floor(selected)) ||
        any(selected < 1L) || any(selected > length(trait_names))) {
      .fd_abort("`trait` contains invalid trait-column indices.")
    }
    selected <- trait_names[as.integer(selected)]
  } else {
    selected <- as.character(selected)
  }
  if (!length(selected) || anyNA(selected) || any(!selected %in% trait_names)) {
    .fd_abort("Requested `trait` is not in the trait table.")
  }

  explicit <- !is.null(ta_axis)
  if (is.null(ta_axis)) {
    if (length(selected) != 1L) {
      .fd_abort(
        "Multivariate trait vectors have no unique total ordering for T-A decomposition. ",
        "Supply `ta_axis` as a trait column or named species score."
      )
    }
    ta_axis <- selected
  }

  species <- rownames(traits)
  if (is.null(species) || anyNA(species) || anyDuplicated(species)) {
    .fd_abort("Trait rows must have unique species names before defining `ta_axis`.")
  }
  if (is.character(ta_axis) && length(ta_axis) == 1L) {
    if (!ta_axis %in% trait_names || !ta_axis %in% names(traits)) {
      .fd_abort("`ta_axis` must name one numeric trait column in the trait table.")
    }
    scores <- traits[[ta_axis]]
    if (!is.numeric(scores) || any(!is.finite(scores))) {
      .fd_abort("The selected `ta_axis` trait column must contain finite numeric scores.")
    }
    return(list(type = "trait_column", name = ta_axis,
                selected_traits = selected, explicit = explicit,
                scores = stats::setNames(as.numeric(scores), species)))
  }
  if (!is.numeric(ta_axis) || is.null(names(ta_axis)) ||
      anyNA(names(ta_axis)) || any(!nzchar(names(ta_axis))) ||
      anyDuplicated(names(ta_axis))) {
    .fd_abort(
      "For multivariate traits, `ta_axis` must be one trait-column name or a named numeric species-score vector."
    )
  }
  if (!all(species %in% names(ta_axis))) {
    missing <- species[!species %in% names(ta_axis)]
    .fd_abort("`ta_axis` is missing scores for species: ", paste(missing, collapse = ", "), ".")
  }
  scores <- as.numeric(ta_axis[species])
  if (any(!is.finite(scores))) {
    .fd_abort("`ta_axis` species scores must be finite and non-missing.")
  }
  list(type = "species_scores", name = NULL,
       selected_traits = selected, explicit = explicit,
       scores = stats::setNames(scores, species))
}

.fd_prepare_pair <- function(reference, focal, x) {
  ref <- .fd_state(x, reference$series, reference$time)
  foc <- .fd_state(x, focal$series, focal$time)
  species <- union(colnames(ref$comm), colnames(foc$comm))
  fallback_traits <- .fd_state_traits(x, species)
  fallback_traits <- fallback_traits[species, , drop = FALSE]
  ref_traits <- fallback_traits
  focal_traits <- fallback_traits
  common_ref <- intersect(rownames(ref$traits), species)
  common_focal <- intersect(rownames(foc$traits), species)
  ref_traits[common_ref, ] <- ref$traits[common_ref, , drop = FALSE]
  focal_traits[common_focal, ] <- foc$traits[common_focal, , drop = FALSE]
  attr(fallback_traits, "species_col") <- x$species_col
  out <- list(ref = .fd_union_state(ref, species, ref_traits),
              focal = .fd_union_state(foc, species, focal_traits),
              traits = fallback_traits, species = species,
              ta_axis = NULL)
  if (length(x$trait_names) == 1L) {
    out$ta_axis <- .fd_resolve_ta_axis(
      ta_axis = NULL, traits = fallback_traits,
      trait_names = x$trait_names, analysed_trait = x$trait_names
    )
  }
  out
}

.fd_subset_state <- function(pair, active, x, pairing_order = NULL,
                               target_order = NULL, axis_order = NULL) {
  if (is.null(pair$ta_axis) || is.null(pair$ta_axis$scores)) {
    .fd_abort(
      "Multivariate trait vectors have no unique total ordering for T-A decomposition. ",
      "Supply `ta_axis` as a trait column or named species score."
    )
  }
  components <- .fd_components_for(x)
  active <- intersect(active, components)
  ref <- pair$ref
  foc <- pair$focal
  if (!length(active)) return(ref)
  if (setequal(active, components)) return(foc)
  use_focal_N <- "N" %in% active
  use_focal_agg <- "aggregation" %in% active
  use_focal_sad <- "SAD" %in% active
  use_focal_traits <- "traits" %in% active
  use_focal_itv <- "ITV" %in% active
  use_focal_ta <- "TA" %in% active
  N <- if (use_focal_N) sum(foc$comm) else sum(ref$comm)
  count_data <- isTRUE(ref$count_data) && isTRUE(foc$count_data)
  if (sum(ref$comm) <= 0 || sum(foc$comm) <= 0) {
    .fd_abort("Temporal decomposition requires positive abundance at both endpoints.")
  }
  sad_source <- if (use_focal_sad) foc else ref
  aggregation_source <- if (use_focal_agg) foc else ref
  std_source <- if (use_focal_traits) foc else ref
  pairing_source <- if (use_focal_ta) foc else ref
  if (use_focal_N && use_focal_agg && use_focal_sad) {
    # These three focal properties uniquely specify the observed focal
    # community structure. Re-running IPF here only approximates a matrix
    # already available exactly and can introduce tiny positive structural
    # zeros. Trait, ITV and TA margins are still constructed below.
    comm <- foc$comm[, pair$species, drop = FALSE]
    construction <- list(
      method = "exact_focal_community",
      row_targets = rowSums(comm),
      col_targets = colSums(comm),
      max_margin_error = 0
    )
  } else {
    target_abundance <- .fd_abundance_slots(
      sad_source, pair$species, N, count_data
    )
    pattern <- .fd_ranked_aggregation_pattern(
      aggregation_source, target_abundance
    )
    # Aggregation profiles, SAD values and the Trait-Abundance pairing are all
    # transported through abundance-rank slots. Species labels therefore do not
    # carry the association that the TA transition is intended to restore.
    comm <- if (count_data) {
      .fd_integer_counterfactual(
        pattern, rowSums(aggregation_source$comm), target_abundance, N
      )
    } else {
      .fd_continuous_counterfactual(
        pattern, rowSums(aggregation_source$comm), target_abundance, N
      )
    }
    construction <- attr(comm, "construction")
    attr(comm, "construction") <- NULL
  }
  traits <- .fd_ranked_trait_margin(
    std_source, pairing_source, colSums(comm),
    pairing_order = pairing_order, target_order = target_order,
    ta_axis_scores = pair$ta_axis$scores, axis_order = axis_order
  )
  itv_source <- if (use_focal_itv) foc else ref
  itv <- .fd_ranked_itv_margin(
    itv_source, traits, colSums(comm), target_order = target_order,
    ta_axis_scores = pair$ta_axis$scores, axis_order = axis_order
  )
  out <- foc
  out$comm <- comm
  # The intermediate community may follow the reference sample layout when
  # aggregation is not yet switched to the focal state. Keep metadata aligned
  # with the rows of the community matrix; this matters when the two dates
  # contain different numbers of sampled sites.
  out$meta <- if (use_focal_agg) foc$meta else ref$meta
  out$traits <- traits
  out$itv <- itv
  out$count_data <- count_data
  out$N <- sum(comm)
  out$species <- colnames(comm)
  out$construction <- c(construction, list(
    SAD_representation = "unlabelled_abundance_rank_slots",
    STD_representation = "unlabelled_complete_trait_rows",
    ITV_representation = "unlabelled_complete_residual-vector distributions conditional on trait-centroid ranks",
    TA_representation = paste0(
      "abundance rank to complete-trait-row permutation along declared scalar axis (",
      if (identical(pair$ta_axis$type, "trait_column")) pair$ta_axis$name else "species scores",
      ")"
    ),
    sources = c(
      N = if (use_focal_N) "focal" else "reference",
      aggregation = if (use_focal_agg) "focal" else "reference",
      SAD = if (use_focal_sad) "focal" else "reference",
      STD = if (use_focal_traits) "focal" else "reference",
      ITV = if (use_focal_itv) "focal" else "reference",
      TA = if (use_focal_ta) "focal" else "reference"
    )
  ))
  out
}

.fd_subset_state_variants <- function(pair, active, x,
                                        selected_components = .fd_components_for(x),
                                        tie_method = c("average", "random", "legacy"),
                                        n_tie_perm = 200L,
                                        max_exact_tie_perm = 10000L,
                                        tie_seed = NULL) {
  tie_method <- match.arg(tie_method)
  components <- .fd_components_for(x)
  if (!length(active) || setequal(intersect(active, components), selected_components) &&
      setequal(selected_components, components)) {
    observed <- if (!length(active)) pair$ref else pair$focal
    return(list(states = list(observed), diagnostics = list(
      method = "observed", tied = FALSE, exact = TRUE, total = 1,
      n_used = 1L, group_sizes = integer(),
      abundance_tie_group_sizes = integer(), trait_axis_tie_group_sizes = integer(),
      exact_or_monte_carlo = "exact", n_admissible_pairings = 1,
      n_evaluated_pairings = 1L, tie_seed = tie_seed
    )))
  }
  base <- .fd_subset_state(pair, active, x)

  pairing_source <- if ("TA" %in% active) pair$focal else pair$ref
  std_source <- if ("traits" %in% active) pair$focal else pair$ref
  source_abundance <- colSums(pairing_source$comm)
  target_abundance <- colSums(base$comm)
  source_key <- .fd_species_tie_key(pairing_source,
    pair$ta_axis$scores[names(source_abundance)], pair$traits)
  target_key <- .fd_species_tie_key(base,
    pair$ta_axis$scores[names(target_abundance)], pair$traits)
  source_ties <- .fd_tie_components(source_abundance, source_key)
  target_ties <- .fd_tie_components(target_abundance, target_key)
  abundance_total <- if (is.infinite(source_ties$count) || is.infinite(target_ties$count) ||
                        source_ties$count > max_exact_tie_perm / target_ties$count) Inf else
    source_ties$count * target_ties$count

  relevant_species <- unique(c(
    names(source_abundance)[source_abundance > 0],
    names(colSums(std_source$comm))[colSums(std_source$comm) > 0],
    names(target_abundance)[target_abundance > 0],
    names(colSums(pair$ref$comm))[colSums(pair$ref$comm) > 0],
    names(colSums(pair$focal$comm))[colSums(pair$focal$comm) > 0]
  ))
  axis_scores <- pair$ta_axis$scores[relevant_species]
  axis_ties <- .fd_axis_tie_components(axis_scores, pair$traits)
  effective_abundance_total <- if (tie_method == "legacy") 1 else abundance_total
  combined_total <- if (is.infinite(effective_abundance_total) || is.infinite(axis_ties$count) ||
                        effective_abundance_total > max_exact_tie_perm / axis_ties$count) Inf else
    effective_abundance_total * axis_ties$count
  exact <- tie_method == "average" && is.finite(combined_total) &&
    combined_total <= max_exact_tie_perm
  axis_method <- if (tie_method == "legacy") "average" else tie_method
  axis_exact <- axis_method == "average" && is.finite(combined_total) &&
    combined_total <= max_exact_tie_perm
  abundance_plan <- .fd_tie_plan(
    source_abundance = source_abundance, target_abundance = target_abundance,
    source_key = source_key, target_key = target_key,
    tie_method = tie_method, n_tie_perm = n_tie_perm,
    max_exact_tie_perm = if (exact) max_exact_tie_perm else 1L,
    seed = tie_seed
  )
  axis_seed <- if (is.null(tie_seed)) NULL else
    as.integer((as.double(tie_seed) + 7919) %% .Machine$integer.max)
  axis_orders <- .fd_axis_tie_orders(axis_ties, axis_exact,
    if (axis_method == "random") 1L else n_tie_perm, axis_seed)
  if (exact) {
    order_pairs <- lapply(abundance_plan$orders, function(ab_order) {
      lapply(axis_orders, function(ax) list(abundance = ab_order, axis = ax))
    })
    order_pairs <- unlist(order_pairs, recursive = FALSE)
  } else {
    n_pairings <- if (tie_method == "random") 1L else as.integer(n_tie_perm)
    order_pairs <- lapply(seq_len(n_pairings), function(i) list(
      abundance = abundance_plan$orders[[((i - 1L) %% length(abundance_plan$orders)) + 1L]],
      axis = axis_orders[[((i - 1L) %% length(axis_orders)) + 1L]]
    ))
  }

  states <- lapply(order_pairs, function(order_pair) {
    st <- base
    pairing_source <- if ("TA" %in% active) pair$focal else pair$ref
    st$traits <- .fd_ranked_trait_margin(
      std_source, pairing_source, colSums(base$comm),
      pairing_order = order_pair$abundance$source,
      target_order = order_pair$abundance$target,
      ta_axis_scores = pair$ta_axis$scores, axis_order = order_pair$axis
    )
    itv_source <- if ("ITV" %in% active) pair$focal else pair$ref
    st$itv <- .fd_ranked_itv_margin(
      itv_source, st$traits, colSums(base$comm),
      target_order = order_pair$abundance$target,
      ta_axis_scores = pair$ta_axis$scores, axis_order = order_pair$axis
    )
    st
  })
  tied <- source_ties$count > 1 || target_ties$count > 1 || axis_ties$tied
  diagnostics <- list(
    method = if (!tied) "none" else if (exact) "average_exact" else "average_monte_carlo",
    tied = tied, exact = exact, total = combined_total,
    n_used = length(order_pairs),
    group_sizes = c(source_ties$group_sizes, target_ties$group_sizes,
                    axis_ties$group_sizes),
    source_group_sizes = source_ties$group_sizes,
    target_group_sizes = target_ties$group_sizes,
    abundance_tie_group_sizes = c(source_ties$group_sizes, target_ties$group_sizes),
    trait_axis_tie_group_sizes = axis_ties$group_sizes[axis_ties$group_sizes > 1L],
    exact_or_monte_carlo = if (exact) "exact" else "monte_carlo",
    n_admissible_pairings = combined_total,
    n_evaluated_pairings = length(order_pairs),
    tie_seed = if (exact || !tied) NULL else tie_seed,
    axis_tied = axis_ties$tied
  )
  list(states = states, diagnostics = diagnostics)
}

.fd_eval_tie_curve <- function(state_info, metric, effort, effort_unit, trait,
                                 n_rep, seed, spat_algo = NULL) {
  states <- state_info$states
  reps <- lapply(states, function(state) {
    # The same random rarefaction/sample order is used for every tie variant.
    # This is a common-random-number coupling and reduces noise in contrasts.
    .fd_eval_curve_state(state, metric, effort, effort_unit, trait, n_rep,
                           seed, spat_algo, TRUE)
  })
  averaged <- Reduce("+", reps) / length(reps)
  variant_means <- do.call(rbind, lapply(reps, colMeans))
  mc_se <- if (isTRUE(state_info$diagnostics$exact) || nrow(variant_means) < 2L) {
    rep(0, ncol(variant_means))
  } else {
    apply(variant_means, 2L, stats::sd) / sqrt(nrow(variant_means))
  }
  list(replicates = averaged, diagnostics = c(
    state_info$diagnostics,
    list(mc_se = mc_se)
  ))
}

.fd_effort_grid <- function(pair, effort = NULL, effort_unit = "full") {
  if (!is.null(effort)) return(as.numeric(effort))
  if (effort_unit == "full") return(0)
  if (effort_unit == "sample") return(seq_len(min(nrow(pair$ref$comm), nrow(pair$focal$comm))))
  seq_len(min(sum(pair$ref$comm), sum(pair$focal$comm)))
}

.fd_eval_scenario <- function(states, metric, effort, effort_unit, trait, n_rep, seed) {
  vapply(seq_along(states), function(i) {
    .fd_eval_state(states[[i]], metric, effort, effort_unit, trait, n_rep,
                     if (is.null(seed)) NULL else seed + i * 1009L)
  }, numeric(1))
}

.fd_curve_quantile <- function(x, probability) {
  apply(x, 2L, function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) NA_real_ else stats::quantile(z, probability, names = FALSE)
  })
}

.fd_permutations <- function(seed = NULL, n_perm = 120L, components = .fd_components) {
  if (length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 1 || n_perm != floor(n_perm)) {
    .fd_abort("`n_perm` must be one positive integer.")
  }
  .fd_with_seed(seed, function() {
    n_orders <- factorial(length(components))
    all <- if (n_perm == n_orders) {
      rec <- function(x) {
        if (!length(x)) return(list(character()))
        unlist(lapply(seq_along(x), function(i) lapply(rec(x[-i]), function(y) c(x[i], y))), recursive = FALSE)
      }
      rec(components)
    } else NULL
    if (!is.null(all)) return(all)
    replicate(as.integer(n_perm), sample(components), simplify = FALSE)
  })
}

#' Compute Shapley path-sensitivity summaries for a decomposition
#'
#' Order variation is structural path sensitivity, not sampling uncertainty
#' and not a confidence interval.
#' @param object An `fundecomp_result` produced by `decompose_fmetric_change()`.
#' @param exact Enumerate every replacement order. Defaults to `TRUE`.
#' @param n_orders Number of randomly sampled orders when `exact = FALSE`.
#' @param seed Optional seed for sampled orders.
#' @return An `fundecomp_result` object with exact or sampled Shapley summaries.
#' @export
shapley_decomposition <- function(object, exact = TRUE, n_orders = NULL, seed = NULL) {
  if (!inherits(object, "fundecomp_result") || is.null(object$.input)) {
    .fd_abort("`object` must be a decomposition result with retained endpoint input.")
  }
  components <- object$components
  total <- factorial(length(components))
  if (isTRUE(exact)) n_orders <- total else {
    if (is.null(n_orders) || length(n_orders) != 1L || !is.finite(n_orders) ||
        n_orders < 1 || n_orders != floor(n_orders)) {
      .fd_abort("Supply a positive integer `n_orders` when `exact = FALSE`.")
    }
  }
  args <- object$.call
  args$method <- "shapley"
  args$n_perm <- as.integer(n_orders)
  args$seed <- if (is.null(seed)) args$seed else seed
  do.call(.fd_decompose_impl, c(list(x = object$.input), args))
}

.fd_decompose_impl <- function(x, reference, focal, metric = "FRic", effort = NULL,
                              effort_unit = "full", method = c("mobr_path", "shapley"),
                              path = NULL, trait = NULL,
                              ta_axis = NULL,
                              n_perm = 1L, seed = NULL, spat_algo = NULL,
                              tie_method = c("average", "random", "legacy"),
                              n_tie_perm = 200L, max_exact_tie_perm = 10000L,
                              tie_seed = NULL, components = NULL) {
  if (!inherits(x, "fundecomp_input")) .fd_abort("`x` must be an `fundecomp_input` object.")
  if (is.null(reference$time) && !is.null(reference$state)) reference$time <- reference$state
  if (is.null(focal$time) && !is.null(focal$state)) focal$time <- focal$state
  if (is.null(reference$time) || is.null(focal$time)) {
    .fd_abort("`reference` and `focal` must identify observed endpoints with `state` (or legacy `time`).")
  }
  method <- match.arg(method)
  all_components <- .fd_components_for(x)
  if (is.null(components)) components <- all_components
  components <- as.character(components)
  components[components == "STD"] <- "traits"
  components[components == "Aggregation"] <- "aggregation"
  if (!length(components) || anyNA(components) || anyDuplicated(components) ||
      any(!components %in% all_components)) {
    .fd_abort("`components` must contain one or more unique supported components: ",
              paste(all_components, collapse = ", "), ".")
  }
  components <- all_components[all_components %in% components]
  if (is.null(path)) path <- components
  path <- as.character(path)
  path[path == "STD"] <- "traits"
  path[path == "Aggregation"] <- "aggregation"
  if (!setequal(path, components) || length(path) != length(components)) {
    .fd_abort("`path` must contain each selected component exactly once: ", paste(components, collapse = ", "), ".")
  }
  tie_method <- match.arg(tie_method)
  if (tie_method == "legacy") {
    warning("tie_method = legacy is label-dependent when abundance ties occur.",
            call. = FALSE)
  }
  if (length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 1 || n_perm != floor(n_perm)) {
    .fd_abort("`n_perm` must be one positive integer.")
  }
  n_perm <- as.integer(n_perm)
  pair <- .fd_prepare_pair(reference, focal, x)
  metric_info <- .fd_get_metric(metric)
  pair$ta_axis <- .fd_resolve_ta_axis(
    ta_axis = ta_axis, traits = pair$traits, trait_names = x$trait_names,
    analysed_trait = trait
  )
  effort <- .fd_effort_grid(pair, effort, effort_unit)
  effective_seed <- .fd_effective_seed(seed)
  tie_effective_seed <- if (is.null(tie_seed)) effective_seed else .fd_effective_seed(tie_seed)
  seed_out <- effective_seed
  if (effort_unit == "individual" && any(effort != floor(effort))) .fd_abort("Individual effort must be integer.")
  if (method == "mobr_path") {
    active <- list(character())
    for (j in seq_along(path)) active[[j + 1L]] <- c(active[[j]], path[j])
    state_infos <- lapply(active, .fd_subset_state_variants, pair = pair, x = x,
                          selected_components = components,
                          tie_method = tie_method, n_tie_perm = n_tie_perm,
                          max_exact_tie_perm = max_exact_tie_perm,
                          tie_seed = tie_effective_seed)
    states <- lapply(state_infos, function(info) info$states[[1L]])
    evaluated <- lapply(state_infos, function(info)
      .fd_eval_tie_curve(info, metric_info, effort, effort_unit, trait,
                           n_perm, effective_seed, spat_algo))
    replicate_values <- lapply(evaluated, function(z) z$replicates)
    tie_diagnostics <- lapply(evaluated, function(z) z$diagnostics)
    vals <- lapply(replicate_values, colMeans)
    reference_metric_values <- vals[[1L]]
    partial_metric_values <- vals[[length(vals)]]
    components <- path
    contrast_replicates <- stats::setNames(
      lapply(seq_along(path), function(j)
        replicate_values[[j + 1L]] - replicate_values[[j]]),
      path
    )
    estimates <- lapply(contrast_replicates, colMeans)
    lower <- lapply(contrast_replicates, .fd_curve_quantile, probability = 0.025)
    upper <- lapply(contrast_replicates, .fd_curve_quantile, probability = 0.975)
    endpoint_replicates <- replicate_values[[length(replicate_values)]] - replicate_values[[1L]]
    endpoint <- colMeans(endpoint_replicates)
    if (setequal(components, all_components)) {
      observed_change <- endpoint
    } else {
      observed_info <- .fd_subset_state_variants(pair, all_components, x,
        selected_components = all_components, tie_method = tie_method,
        n_tie_perm = n_tie_perm, max_exact_tie_perm = max_exact_tie_perm,
        tie_seed = tie_effective_seed)
      observed_replicates <- .fd_eval_tie_curve(observed_info, metric_info, effort,
        effort_unit, trait, n_perm, effective_seed, spat_algo)$replicates
      observed_change <- colMeans(observed_replicates - replicate_values[[1L]])
    }
    recon <- Reduce("+", estimates) - endpoint
    rows <- lapply(seq_along(components), function(j) data.frame(
      series = as.character(reference$series), time_reference = as.character(reference$time),
      time_focal = as.character(focal$time), metric = metric_info$name,
      effort = effort, component = components[j],
      component_label = unname(.fd_component_labels[components[j]]),
      estimate = estimates[[j]],
      lower = lower[[j]], upper = upper[[j]], method = method,
      n_perm = n_perm, seed = seed_out,
      reconstruction_error = recon, stringsAsFactors = FALSE
    ))
    curves <- do.call(rbind, rows)
    scenario_labels <- c("reference", paste0("after_", path))
    scenario_curves <- do.call(rbind, lapply(seq_along(vals), function(i) {
      data.frame(state = scenario_labels[i], state_index = i - 1L,
                 effort = effort, value = vals[[i]],
                 stringsAsFactors = FALSE)
    }))
    permutations <- NULL
  } else {
    perms <- .fd_permutations(effective_seed, max(1L, as.integer(n_perm)), components)
    contrib <- array(NA_real_, dim = c(length(perms), length(components), length(effort)),
                     dimnames = list(NULL, components, NULL))
    coalition_key <- function(active) {
      active <- components[components %in% active]
      if (!length(active)) "empty" else paste(active, collapse = "+")
    }
    # Every unique coalition is evaluated once and reused across replacement
    # orders; this is exactly equivalent to repeated evaluation because every
    # stochastic state receives the same effective seed.
    coalition_sets <- list(character())
    for (permutation in perms) {
      for (j in seq_along(permutation)) {
        coalition_sets[[length(coalition_sets) + 1L]] <- permutation[seq_len(j)]
      }
    }
    coalition_keys <- vapply(coalition_sets, coalition_key, character(1))
    keep_coalition <- !duplicated(coalition_keys)
    coalition_sets <- coalition_sets[keep_coalition]
    coalition_keys <- coalition_keys[keep_coalition]
    coalition_infos <- lapply(coalition_sets, function(active) {
      info <- .fd_subset_state_variants(
        pair, active, x, selected_components = components, tie_method = tie_method, n_tie_perm = n_tie_perm,
        max_exact_tie_perm = max_exact_tie_perm, tie_seed = tie_effective_seed
      )
      evaluated <- .fd_eval_tie_curve(
        info, metric_info, effort, effort_unit, trait, 1L, effective_seed,
        spat_algo
      )
      list(value = as.numeric(evaluated$replicates),
           diagnostics = evaluated$diagnostics)
    })
    coalition_values <- lapply(coalition_infos, function(z) z$value)
    names(coalition_values) <- coalition_keys
    names(coalition_infos) <- coalition_keys
    for (p in seq_along(perms)) {
      active <- character()
      prev <- coalition_values[[coalition_key(active)]]
      for (j in seq_along(perms[[p]])) {
        active <- c(active, perms[[p]][j])
        now <- coalition_values[[coalition_key(active)]]
        contrib[p, perms[[p]][j], ] <- now - prev
        prev <- now
      }
    }
    estimates <- apply(contrib, c(2, 3), mean, na.rm = TRUE)
    lower <- apply(contrib, c(2, 3), stats::quantile, probs = 0.025, na.rm = TRUE)
    upper <- apply(contrib, c(2, 3), stats::quantile, probs = 0.975, na.rm = TRUE)
    endpoint <- coalition_values[[coalition_key(components)]] -
      coalition_values[[coalition_key(character())]]
    reference_metric_values <- coalition_values[[coalition_key(character())]]
    partial_metric_values <- coalition_values[[coalition_key(components)]]
    if (setequal(components, all_components)) {
      observed_change <- endpoint
    } else {
      observed_info <- .fd_subset_state_variants(pair, all_components, x,
        selected_components = all_components, tie_method = tie_method,
        n_tie_perm = n_tie_perm, max_exact_tie_perm = max_exact_tie_perm,
        tie_seed = tie_effective_seed)
      observed_value <- as.numeric(.fd_eval_tie_curve(observed_info, metric_info,
        effort, effort_unit, trait, 1L, effective_seed, spat_algo)$replicates)
      observed_change <- observed_value - coalition_values[[coalition_key(character())]]
    }
    recon <- colSums(estimates) - endpoint
    curves <- do.call(rbind, lapply(seq_along(components), function(j) data.frame(
      series = as.character(reference$series), time_reference = as.character(reference$time),
      time_focal = as.character(focal$time), metric = metric_info$name,
      effort = effort, component = components[j],
      component_label = unname(.fd_component_labels[components[j]]),
      estimate = estimates[j, ],
      lower = lower[j, ], upper = upper[j, ], method = method,
      n_perm = n_perm, seed = seed_out,
      reconstruction_error = recon, stringsAsFactors = FALSE
    )))
    states <- NULL
    scenario_curves <- NULL
    contrast_replicates <- NULL
    endpoint_replicates <- NULL
    permutations <- perms
    tie_diagnostics <- coalition_infos
  }
  curves$interval_type <- if (method == "shapley") {
    "structural_path_sensitivity"
  } else {
    "paired_monte_carlo_contrast_quantiles"
  }
  curves$selected_components <- I(rep(list(components), nrow(curves)))
  curves$reference_metric <- rep(reference_metric_values, times = length(components))
  curves$partial_counterfactual_metric <- rep(partial_metric_values, times = length(components))
  curves$observed_focal_metric <- rep(reference_metric_values + observed_change, times = length(components))
  curves$partial_change <- rep(endpoint, times = length(components))
  curves$observed_total_change <- rep(observed_change, times = length(components))
  curves$remaining_difference <- rep(observed_change - endpoint, times = length(components))
  curves$reconstruction_residual <- rep(recon, times = length(components))
  out <- list(curves = curves, states = states,
              states_are_representatives = !is.null(states),
              reference = reference, focal = focal,
              scenario_curves = scenario_curves,
              metric = metric_info$name, effort = effort, effort_unit = effort_unit,
              method = method, path = path,
              components = components,
              full_decomposition = setequal(components, all_components),
              pool_semantics = list(
                sad = "descriptive_only; observed endpoint SADs define direct decomposition"
              ),
              ta_axis = pair$ta_axis,
              constraints = NULL,
              permutations = permutations,
              contrast_replicates = contrast_replicates,
              endpoint_replicates = endpoint_replicates,
              diagnostics = list(
                tie = list(
                  method = tie_method, n_tie_perm = as.integer(n_tie_perm),
                  max_exact_tie_perm = as.integer(max_exact_tie_perm),
                  seed = tie_effective_seed,
                  any_ties = any(vapply(
                    tie_diagnostics, function(z) isTRUE(z$tied), logical(1)
                  )),
                  states = tie_diagnostics
                ),
                order_sensitivity = if (method == "shapley") list(
                  interpretation = "structural_path_sensitivity_not_sampling_uncertainty",
                  n_orders = length(permutations), exact = length(permutations) == factorial(length(components))
                ) else NULL
              ),
              endpoint = list(reference = pair$ref, focal = pair$focal),
              settings = list(path = path, metric = metric_info$name,
                trait = trait, ta_axis = ta_axis, species_support = pair$species,
                effort_unit = effort_unit, n_rep = n_perm, seed = seed_out),
              .input = x, .call = list(reference = reference, focal = focal,
                metric = metric, effort = effort, effort_unit = effort_unit,
                path = path, trait = trait, ta_axis = ta_axis, seed = seed_out,
                spat_algo = spat_algo, tie_method = tie_method,
                n_tie_perm = as.integer(n_tie_perm),
                max_exact_tie_perm = as.integer(max_exact_tie_perm),
                tie_seed = tie_seed, components = components))
  class(out) <- "fundecomp_result"
  out$constraints <- .fd_transition_constraint_table(out)
  out
}

#' Decompose change between two observed assemblage states
#'
#' The reference and focal states may be times, treatments, habitats, or other
#' ecological contrasts. Hybrid states use endpoint properties only.
#' With `components`, only the selected properties are changed to their focal
#' values; other properties remain at reference. The reported contributions
#' then reconstruct this conditional contrast, which can differ from the
#' observed focal-reference change. For Shapley attribution, selected
#' components are the players in a restricted coalition game.
#' Aggregation denotes sample-level allocation conditional on total abundance
#' and the global SAD; sample totals may change at this transition. Tied
#' abundance ranks and tied values on the declared TA axis are averaged over
#' admissible pairings. See [validate_decomposition()] for machine-checkable
#' reconstruction and margin checks.
#' @param x An `fundecomp_input` object.
#' @param reference Reference endpoint, with `series` and `time` (or `state`).
#' @param focal Focal endpoint, with `series` and `time` (or `state`).
#' @param metric Metric name or custom metric function.
#' @param trait Trait column(s) required by the metric.
#' @param ta_axis Explicit scalar axis for multivariate Trait-Abundance ordering.
#' @param series Reserved for selecting generic series in the input.
#' @param effort Shared effort grid.
#' @param effort_unit Effort scale: `individual`, `sample`, or `full`.
#' @param accumulation Sample accumulation convention.
#' @param method Attribution method: a declared sequential path or exact
#'   Shapley averaging over component orders.
#' @param path Component replacement order.
#' @param components Components to change from reference to focal. Unselected
#'   components remain at reference. The resulting contributions reconstruct
#'   this conditional contrast, not necessarily the observed endpoint change.
#' @param species_support Species alignment; defaults to endpoint union.
#' @param n_rep Number of paired Monte Carlo evaluations.
#' @param seed Optional reproducibility seed.
#' @param tie_method How abundance and TA-axis ties are handled: exact or
#'   Monte Carlo averaging (`"average"`), random resolution, or label-dependent
#'   legacy ordering.
#' @param n_tie_perm Number of pairings sampled when exact tie enumeration is
#'   too large, or when random tie resolution is requested.
#' @param max_exact_tie_perm Maximum combined number of pairings to enumerate.
#' @param tie_seed Optional seed for sampled tie pairings.
#' @param ... Deprecated compatibility arguments. A non-`NULL` `pool` is
#'   rejected because regional-pool baselines are not implemented; pool objects
#'   currently provide descriptive summaries only.
#' @return An `fundecomp_result` decomposition object.
#' @export
decompose_fmetric_change <- function(x, reference, focal, metric = "FRic",
                                    trait = NULL, ta_axis = NULL, series = NULL,
                                    effort = NULL,
                                    effort_unit = c("sample", "individual", "full"),
                                    accumulation = c("nonspatial", "spatial"),
                                    path = NULL,
                                    species_support = "pairwise_union",
                                    method = c("mobr_path", "shapley"),
                                    n_rep = 99L, seed = NULL,
                                    tie_method = c("average", "random", "legacy"),
                                    n_tie_perm = 200L,
                                    max_exact_tie_perm = 10000L,
                                    tie_seed = NULL, components = NULL, ...) {
  effort_unit <- match.arg(effort_unit)
  accumulation <- match.arg(accumulation)
  method <- match.arg(method)
  dots <- list(...)
  if (!is.null(dots$pool)) {
    .fd_abort("`pool` does not change decompositions. Regional-pool SAD/STD baselines are not implemented.")
  }
  deprecated <- intersect(names(dots), c("n_perm", "spat_algo"))
  if (length(deprecated)) warning("Deprecated arguments: ", paste(deprecated, collapse = ", "), call. = FALSE)
  if (!is.null(series)) warning("`series` selection is currently expressed in `reference` and `focal`.", call. = FALSE)
  if (!identical(species_support, "pairwise_union")) {
    .fd_abort("Explicit `species_support` is not yet supported; use the pairwise endpoint union.")
  }
  if (!is.null(dots$n_perm)) n_rep <- dots$n_perm
  if (accumulation == "spatial" && is.null(x$coord_names)) {
    .fd_abort("Spatial accumulation requires coordinates in the input object.")
  }
  spat_algo <- if (accumulation == "spatial") "kNN" else NULL
  if (!is.null(dots$spat_algo)) spat_algo <- dots$spat_algo
  args <- list(x = x, reference = reference, focal = focal, metric = metric,
               effort = effort, effort_unit = effort_unit, method = method,
               path = path, trait = trait, ta_axis = ta_axis,
               n_perm = n_rep, seed = seed, spat_algo = spat_algo,
               tie_method = tie_method, n_tie_perm = n_tie_perm,
               max_exact_tie_perm = max_exact_tie_perm, tie_seed = tie_seed)
  args$components <- components
  do.call(.fd_decompose_impl, args)
}

#' Deprecated v0.1 curve interface
#' @param ... Legacy arguments passed to the core implementation.
#' @export
get_fdelta_curves <- function(...) {
  args <- list(...)
  if (!is.null(args$pool)) {
    .fd_abort("`pool` does not change decompositions. Regional-pool SAD/STD baselines are not implemented.")
  }
  args$pool <- NULL
  do.call(.fd_decompose_impl, args)
}

#' Summarize decomposition curves by component and effort
#' @param x An `fundecomp_result` object.
#' @param effort Optional effort value or values to retain. By default all
#'   efforts are returned; values are never summed across different effort
#'   scales.
#' @param ... Unused arguments.
#' @export
get_fdelta_stats <- function(x, effort = NULL, ...) {
  if (!inherits(x, "fundecomp_result")) .fd_abort("`x` must be an `fundecomp_result` object.")
  z <- x$curves
  if (!is.null(effort)) {
    if (any(!is.finite(effort)) || any(!effort %in% z$effort)) {
      .fd_abort("Requested `effort` is absent from the decomposition curves.")
    }
    z <- z[z$effort %in% effort, , drop = FALSE]
  }
  out <- z[c("component", "component_label", "effort", "estimate", "lower",
             "upper", "reconstruction_error")]
  rownames(out) <- NULL
  out
}
