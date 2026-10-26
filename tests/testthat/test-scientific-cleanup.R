test_that("Aggregation includes sample-total profile changes at fixed N and SAD", {
  comm <- rbind(
    ref_1 = c(sp1 = 1L, sp2 = 0L), ref_2 = c(sp1 = 1L, sp2 = 0L),
    ref_3 = c(sp1 = 0L, sp2 = 2L), foc_1 = c(sp1 = 1L, sp2 = 1L),
    foc_2 = c(sp1 = 1L, sp2 = 0L), foc_3 = c(sp1 = 0L, sp2 = 1L)
  )
  expect_equal(sum(comm[1:3, ]), sum(comm[4:6, ]))
  expect_equal(colSums(comm[1:3, ]), colSums(comm[4:6, ]))
  expect_false(identical(unname(rowSums(comm[1:3, ])), unname(rowSums(comm[4:6, ]))))
  x <- make_fundecomp_input(comm = comm,
                    traits = data.frame(species = c("sp1", "sp2"), mass = c(1, 5)),
                    series = rep("s", 6), time = rep(c(1, 2), each = 3))
  out <- decompose_fmetric_change(
    x, list(series = "s", time = 1), list(series = "s", time = 2),
    metric = "richness", effort = 1, effort_unit = "sample", n_rep = 100,
    seed = 9
  )
  estimate <- stats::setNames(out$curves$estimate, out$curves$component)
  expect_equal(unname(estimate["N"]), 0, tolerance = 1e-12)
  expect_equal(unname(estimate["SAD"]), 0, tolerance = 1e-12)
  expect_gt(abs(unname(estimate["aggregation"])), 0)
  expect_lt(max(abs(out$curves$reconstruction_error)), 1e-12)
})

test_that("transition constraints are explicit and validator detects corruption", {
  comm <- rbind(r1 = c(a = 8L, b = 2L, c = 1L), r2 = c(a = 1L, b = 3L, c = 7L),
                f1 = c(a = 2L, b = 6L, c = 3L), f2 = c(a = 5L, b = 2L, c = 4L))
  x <- make_fundecomp_input(comm = comm,
                    traits = data.frame(species = c("a", "b", "c"), mass = c(1, 4, 8)),
                    series = rep("s", 4), time = rep(c(1, 2), each = 2))
  out <- decompose_fmetric_change(
    x, list(series = "s", time = 1), list(series = "s", time = 2),
    metric = "trait_sum", trait = "mass", effort_unit = "full", seed = 11
  )
  expect_true(all(c("transition", "property", "action", "checked", "passed") %in%
                    names(out$constraints)))
  expect_true(any(out$constraints$action == "conditioned_upon"))
  validation <- validate_decomposition(out)
  expect_true(validation$passed)
  expect_true(all(validation$check$passed[!is.na(validation$check$passed)]))

  corrupted <- out
  corrupted$states[[2L]]$comm[1, 1] <- corrupted$states[[2L]]$comm[1, 1] + 1
  expect_false(validate_decomposition(corrupted)$passed)
})

make_validator_fixture <- function(ref, foc, traits, count_data = TRUE, itv = NULL) {
  comm <- rbind(ref, foc)
  nr <- nrow(ref); nf <- nrow(foc)
  x <- make_fundecomp_input(comm = comm, traits = traits,
    series = rep("s", nr + nf), state = rep(c("ref", "foc"), c(nr, nf)),
    count_data = count_data, trait_observations = itv)
  decompose_fmetric_change(x, list(series = "s", time = "ref"),
    list(series = "s", time = "foc"), metric = "CWM",
    trait = names(traits)[2L], effort_unit = "full", seed = 23)
}

test_that("N transition checks deterministically scaled unlabelled SAD and rows", {
  ref <- rbind(r1 = c(a = 5, b = 3, c = 2), r2 = c(a = 5, b = 2, c = 3))
  foc <- rbind(f1 = c(a = 10, b = 5, c = 5), f2 = c(a = 10, b = 5, c = 5))
  traits <- data.frame(species = c("a", "b", "c"), value = c(1, 3, 7))
  out <- make_validator_fixture(ref, foc, traits)
  v <- validate_decomposition(out)
  expect_true(v$check$passed[match("N_SAD", v$check$check)])
  expect_true(v$check$passed[match("N_sample_allocation", v$check$check)])
  expect_equal(sum(out$states[[2L]]$comm), 40)

  corrupted <- out
  corrupted$states[[2L]]$comm[1, "a"] <- corrupted$states[[2L]]$comm[1, "a"] + 1L
  corrupted$states[[2L]]$comm[1, "b"] <- corrupted$states[[2L]]$comm[1, "b"] - 1L
  expect_false(validate_decomposition(corrupted)$check$passed[
    match("N_SAD", validate_decomposition(corrupted)$check$check)])

  bad_rows <- out
  bad_rows$states[[2L]]$comm[1, "a"] <- bad_rows$states[[2L]]$comm[1, "a"] + 1L
  bad_rows$states[[2L]]$comm[2, "a"] <- bad_rows$states[[2L]]$comm[2, "a"] - 1L
  expect_false(validate_decomposition(bad_rows)$check$passed[
    match("N_sample_allocation", validate_decomposition(bad_rows)$check$check)])
})

test_that("N transition validates continuous row and SAD scaling", {
  ref <- rbind(r1 = c(a = 5, b = 3, c = 2), r2 = c(a = 5, b = 2, c = 3))
  foc <- rbind(f1 = c(a = 10, b = 5, c = 5), f2 = c(a = 10, b = 5, c = 5))
  traits <- data.frame(species = c("a", "b", "c"), value = c(1, 3, 7))
  out <- make_validator_fixture(ref, foc, traits, count_data = FALSE)
  checks <- validate_decomposition(out)$check
  expect_true(checks$passed[match("N_SAD", checks$check)])
  expect_true(checks$passed[match("N_sample_allocation", checks$check)])
})

test_that("TA community validation is invariant to labels and tied column order", {
  unique_ref <- rbind(ur1 = c(a = 6, b = 3, c = 1), ur2 = c(a = 0, b = 0, c = 0))
  unique_foc <- rbind(uf1 = c(a = 4, b = 3, c = 2), uf2 = c(a = 0, b = 0, c = 0))
  unique_out <- make_validator_fixture(unique_ref, unique_foc,
    data.frame(species = c("a", "b", "c"), value = c(1, 5, 8)))
  unique_checks <- validate_decomposition(unique_out)$check
  expect_true(unique_checks$passed[match("TA_community_fixed", unique_checks$check)])
  ta_constraint <- unique_out$constraints[
    unique_out$constraints$transition == "TA" &
      unique_out$constraints$property == "sample-level allocation", , drop = FALSE]
  expect_true(nrow(ta_constraint) == 1L && ta_constraint$checked && ta_constraint$passed)

  ref <- rbind(r1 = c(a = 4, b = 1, c = 1), r2 = c(a = 0, b = 1, c = 1))
  foc <- rbind(f1 = c(a = 2, b = 1, c = 1), f2 = c(a = 2, b = 1, c = 1))
  traits <- data.frame(species = c("a", "b", "c"), value = c(1, 5, 5))
  out <- make_validator_fixture(ref, foc, traits)
  checks <- validate_decomposition(out)$check
  expect_true(checks$passed[match("TA_community_fixed", checks$check)])

  # Counterfactual scaling can leave tiny numerical dust in the zero cells of
  # one tied representative. It must not split the abundance class or change
  # the canonical ordering of its sample profiles.
  dust <- out
  ta_before <- length(dust$states) - 1L
  tied <- which(abs(colSums(dust$states[[ta_before]]$comm) - 2) < 1e-8)
  if (length(tied) >= 2L) {
    z <- dust$states[[ta_before]]$comm
    tied_profiles <- z[, tied, drop = FALSE]
    tied_profiles[tied_profiles == 0] <- 1e-12
    z[, tied] <- tied_profiles
    dust$states[[ta_before]]$comm <- z
    dust_check <- validate_decomposition(dust)$check
    expect_true(dust_check$passed[match("TA_community_fixed", dust_check$check)])
  }

  renamed <- out
  for (i in seq_along(renamed$states)) {
    colnames(renamed$states[[i]]$comm) <- paste0("renamed_", seq_len(ncol(renamed$states[[i]]$comm)))
  }
  expect_true(validate_decomposition(renamed)$check$passed[
    match("TA_community_fixed", validate_decomposition(renamed)$check$check)])

  reordered <- out
  for (i in seq_along(reordered$states))
    reordered$states[[i]]$comm <- reordered$states[[i]]$comm[, rev(seq_len(ncol(reordered$states[[i]]$comm))), drop = FALSE]
  expect_true(validate_decomposition(reordered)$check$passed[
    match("TA_community_fixed", validate_decomposition(reordered)$check$check)])

  # Alter a 2x2 allocation cycle at the TA transition while retaining both
  # margins. The rank-space sample profiles must detect this change.
  corrupted <- out
  ta_before <- length(out$path)
  j <- ta_before
  m <- corrupted$states[[j]]$comm
  tied <- which(colSums(m) == 2)
  if (nrow(m) >= 2L && length(tied) >= 2L && m[1, tied[1L]] > 0 && m[2, tied[2L]] > 0) {
    m[1, tied[1L]] <- m[1, tied[1L]] - 1
    m[2, tied[2L]] <- m[2, tied[2L]] - 1
    m[1, tied[2L]] <- m[1, tied[2L]] + 1
    m[2, tied[1L]] <- m[2, tied[1L]] + 1
    corrupted$states[[j]]$comm <- m
    cv <- validate_decomposition(corrupted)
    expect_false(cv$check$passed[match("TA_community_fixed", cv$check$check)])
  }
})

test_that("SAD to ITV transport is marked approximate across richness change", {
  # The constraint table and validator must not claim a direct ITV identity
  # where systematic transport maps between different positive richness.
  ref <- rbind(r1 = c(a = 4, b = 2, c = 0), r2 = c(a = 2, b = 0, c = 0))
  foc <- rbind(f1 = c(a = 2, b = 2, c = 2), f2 = c(a = 2, b = 0, c = 0))
  traits <- data.frame(species = c("a", "b", "c"), value = c(1, 3, 7))
  obs <- data.frame(sample = c("r1", "r1", "r2", "f1", "f1", "f2"),
    species = c("a", "b", "a", "a", "c", "a"), value = c(1, 2, 1, 1, 2, 1))
  x <- make_fundecomp_input(comm = rbind(ref, foc), traits = traits,
    series = rep("s", 4), state = rep(c("ref", "foc"), each = 2),
    trait_observations = obs)
  out <- decompose_fmetric_change(x, list(series = "s", time = "ref"),
    list(series = "s", time = "foc"), metric = "FDis", trait = "value",
    effort_unit = "full", seed = 91)
  checks <- validate_decomposition(out)$check
  sad_itv <- checks[checks$check == "SAD_ITV", , drop = FALSE]
  if (nrow(sad_itv)) {
    expect_true(is.na(sad_itv$passed) || sad_itv$status == "approximate_systematic_transport")
  }

  same_ref <- ref; same_ref[, "c"] <- c(1, 1); same_ref[, "a"] <- same_ref[, "a"] - c(1, 1)
  same_comm <- rbind(same_ref, foc); rownames(same_comm) <- paste0("same", seq_len(nrow(same_comm)))
  same_rich <- make_fundecomp_input(comm = same_comm, traits = traits,
    series = rep("s", 4), state = rep(c("ref", "foc"), each = 2),
    trait_observations = obs)
  same_out <- decompose_fmetric_change(same_rich, list(series = "s", time = "ref"),
    list(series = "s", time = "foc"), metric = "FDis", trait = "value",
    effort_unit = "full", seed = 92)
  same_check <- validate_decomposition(same_out)$check
  same_sad_itv <- same_check[same_check$check == "SAD_ITV", , drop = FALSE]
  if (nrow(same_sad_itv)) expect_true(same_sad_itv$passed)

  invalid <- same_out
  sad_state <- which(invalid$path == "SAD")
  pos <- sad_state + 1L
  species <- names(invalid$states[[pos]]$itv)[1L]
  invalid$states[[pos]]$itv[[species]][1L, 1L] <-
    invalid$states[[pos]]$itv[[species]][1L, 1L] + 1
  invalid_check <- validate_decomposition(invalid)$check
  sad_itv_row <- invalid_check[invalid_check$check == "SAD_ITV", , drop = FALSE]
  if (nrow(sad_itv_row)) expect_false(sad_itv_row$passed)
})
