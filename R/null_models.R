#' Generate a fixed-margin random-aggregation null model
#'
#' The randomizer preserves total abundance, every sample total and every
#' species total, while randomizing the sample-species allocation. The
#' historical `mobr::get_null_comm(..., "rand_agg")` remains available through
#' `mobr` but is not used as the definition of this fundecomp null.
#'
#' @param comm Community matrix.
#' @param seed Optional seed. The caller's RNG state is restored on exit.
#' @param n_swap Number of attempted two-by-two swaps. The default scales with
#'   the matrix size and is capped to avoid excessive work.
#' @param method Null implementation. The default is the fixed-margin fundecomp
#'   randomizer; `"mobr"` delegates to the historical mobr null for comparison.
#' @return A randomized community matrix with a `constraints` attribute.
#' @export
null_aggregation <- function(comm, seed = NULL, n_swap = NULL,
                                  method = c("fixed_margin", "mobr")) {
  comm <- .fd_as_matrix(comm)
  .fd_check_abundance(comm, TRUE)
  method <- match.arg(method)
  if (method == "mobr") {
    out <- .fd_with_seed(seed, function() mobr::get_null_comm(comm, "rand_agg"))
    attr(out, "constraints") <- list(
      total_abundance = isTRUE(all.equal(sum(comm), sum(out))),
      sample_totals = isTRUE(all.equal(rowSums(comm), rowSums(out))),
      species_totals = isTRUE(all.equal(colSums(comm), colSums(out)))
    )
    return(out)
  }
  if (nrow(comm) < 2L || ncol(comm) < 2L || sum(comm) == 0) return(comm)
  if (is.null(n_swap)) n_swap <- min(1000000L, max(100L, 20L * nrow(comm) * ncol(comm)))
  if (length(n_swap) != 1L || !is.finite(n_swap) || n_swap < 1 || n_swap != floor(n_swap)) {
    .fd_abort("`n_swap` must be one positive integer.")
  }
  out <- .fd_with_seed(seed, function() {
    z <- comm
    for (i in seq_len(as.integer(n_swap))) {
      rr <- sample.int(nrow(z), 2L)
      cc <- sample.int(ncol(z), 2L)
      if (z[rr[1L], cc[1L]] > 0 && z[rr[2L], cc[2L]] > 0) {
        z[rr[1L], cc[1L]] <- z[rr[1L], cc[1L]] - 1L
        z[rr[2L], cc[2L]] <- z[rr[2L], cc[2L]] - 1L
        z[rr[1L], cc[2L]] <- z[rr[1L], cc[2L]] + 1L
        z[rr[2L], cc[1L]] <- z[rr[2L], cc[1L]] + 1L
      } else if (z[rr[1L], cc[2L]] > 0 && z[rr[2L], cc[1L]] > 0) {
        z[rr[1L], cc[2L]] <- z[rr[1L], cc[2L]] - 1L
        z[rr[2L], cc[1L]] <- z[rr[2L], cc[1L]] - 1L
        z[rr[1L], cc[1L]] <- z[rr[1L], cc[1L]] + 1L
        z[rr[2L], cc[2L]] <- z[rr[2L], cc[2L]] + 1L
      }
    }
    z
  })
  attr(out, "constraints") <- list(
    total_abundance = isTRUE(all.equal(sum(comm), sum(out))),
    sample_totals = isTRUE(all.equal(rowSums(comm), rowSums(out))),
    species_totals = isTRUE(all.equal(colSums(comm), colSums(out))),
    spatial_allocation_randomized = !isTRUE(all.equal(comm, out))
  )
  out
}

#' Break the trait-abundance link while preserving both margins
#'
#' Trait dimensions are permuted jointly among species. Therefore the
#' multivariate species-trait distribution is preserved exactly, as are the
#' community abundances; only the species-trait-abundance pairing is changed.
#'
#' @param state A state list containing `comm` and a species-row `traits` table.
#' @param seed Optional seed. The caller's RNG state is restored on exit.
#' @param permutation Optional integer permutation of the species rows, useful
#'   for deterministic tests.
#' @return A list with the randomized `state`, the applied permutation, and
#'   Boolean constraint diagnostics.
#' @export
break_trait_abundance_link <- function(state, seed = NULL, permutation = NULL) {
  if (!is.list(state) || is.null(state$comm) || is.null(state$traits)) {
    .fd_abort("`state` must contain `comm` and `traits`.")
  }
  species <- colnames(state$comm)
  if (is.null(species) || !all(species %in% rownames(state$traits))) {
    .fd_abort("Trait rows must be named by every community species.")
  }
  traits <- state$traits[species, , drop = FALSE]
  n <- nrow(traits)
  if (is.null(permutation)) permutation <- .fd_with_seed(seed, function()
    sample.int(n, n, replace = FALSE)) else permutation <- as.integer(permutation)
  if (length(permutation) != n || !identical(sort(as.integer(permutation)), seq_len(n))) {
    .fd_abort("`permutation` must be a permutation of the trait rows.")
  }
  out <- state
  out$traits <- traits[permutation, , drop = FALSE]
  rownames(out$traits) <- species
  same_margin <- all(vapply(seq_len(ncol(traits)), function(j)
    isTRUE(all.equal(sort(traits[[j]]), sort(out$traits[[j]]))), logical(1)))
  structure(list(
    state = out,
    permutation = as.integer(permutation),
    constraints = list(
      SAD = isTRUE(all.equal(colSums(state$comm), colSums(out$comm))),
      STD = same_margin,
      multivariate_trait_rows = same_margin,
      trait_abundance_pairing = !identical(permutation, seq_len(n))
    )
  ), class = "fundecomp_null_ta")
}
