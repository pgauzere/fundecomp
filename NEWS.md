# fundecomp 5.0.1

* Declared the GitHub source for the `mobr` dependency so development tools and
  GitHub Actions can resolve and install it when it is not available from the
  configured package repositories.

# fundecomp 0.4.0

- Added `components` to endpoint and time-series decomposition interfaces.
  Selected properties change from reference to focal while all other
  properties remain at reference, producing a conditional partial contrast.
- Results now report reference, selected-counterfactual, and observed focal
  metrics separately, with partial change, observed total change, remaining
  difference, and numerical reconstruction residual.
- Defined partial Shapley attribution as a restricted coalition game with
  selected components as players and all other components held at reference.
- Updated decomposition validation and the methodological vignette for
  conditional partial contrasts.
- Verified deterministic full-path estimates against the installed 0.3.0
  package for CWM, FDis, and richness; all component estimates and reconstruction
  errors were identical in the tested examples.

# fundecomp 0.3.0

- Removed the duplicated unreachable branch in `.fd_subset_state()`; a focal
  community specified by focal N, aggregation, and SAD is still reused exactly.
- Removed the unused pool argument from internal state preparation and
  decomposition. Public decomposition interfaces now reject a non-`NULL` pool
  with a clear message. Regional-pool SAD/STD baselines remain unimplemented.
- Made the supported tie-handling controls explicit arguments of
  `decompose_fmetric_change()` and retained them when creating Shapley results.
- Added the methodological vignette and revised the introductory vignette and
  README to distinguish current methods from planned extensions and to document
  the actual metric registry.
- Updated the Martins study script to use the current `accumulation` and
  `n_rep` arguments. Removed its misleading leave-one-out pool sensitivity
  calculation, which did not alter decompositions.
- Excluded rendered vignette HTML files from source package builds; the `.Rmd`
  documents are the maintained sources.

# fundecomp 0.2.0 (development)

- Made `decompose_fmetric_change()` the main endpoint comparison API and added generic `state` labels alongside legacy `time` input.
- Removed pool fields from new decomposition results. `make_fundecomp_pool()` remains temporarily as a deprecated descriptive compatibility helper and is not called by the new decomposition core.
- Reduced the built-in functional registry to CWM, CWV, FDis, FRic, and `trait_sum`; legacy `biomass` and `variance` names warn and map to `trait_sum` and CWV.
- Added `shapley_decomposition()` as a secondary order-sensitivity analysis. Exact defaults enumerate all 120 or 720 orders.
- Made ambiguous multi-state rarefaction selection an error.
- Defined Aggregation as full sample-level allocation conditional on N and the global SAD, including changes in sample-total heterogeneity. The existing numerical transition already followed this definition.
- Added symmetric TA-axis tie averaging, exact through 10,000 combined abundance/trait pairing orders and seeded Monte Carlo above that threshold.
- Added explicit transition constraints and `validate_decomposition()` checks for reconstruction, observed endpoints, integer states, effort grids, N, SAD, and sample-total margins where applicable.
- Documented unchecked properties explicitly rather than describing them as programmatically validated.
- Updated the five standardized case-study run scripts to use the endpoint API and `trait_sum` for mass totals. Full case studies were not rerun.
# fundecomp 5.0.0

* Added `extract_component_changes()` and `plot_component_changes()` as optional
  tools to describe observed ecological component changes independently of
  functional diversity decomposition contributions.
* Component extraction retains rank-based SAD and trait-abundance structure,
  sample-level allocation, unweighted species traits, and optional ITV data.
* Improved component plots with Initial/Final state labels, grouped N and SAD
  bars, ranked within-species spatial allocation profiles, and trait values
  displayed against abundance rank for TA. The combined plot now uses optional
  `patchwork` and collects repeated legends into one shared legend.
* Added numeric rank axes for SAD curves, an optional log rank axis, and an
  optional species-by-site heatmap for observed spatial allocation. TA points
  retain each observed species value while placing abundance ties at shared
  midranks. Aggregation extraction now uses sample IDs from endpoint metadata.
* Long SAD plots now default to two curves on linear rank and logarithmic
  relative-abundance axes. Added a ranked-site spatial heatmap that retains
  each species-specific allocation profile, including zero-abundance sites,
  while discarding site identity. Endpoint colours distinguish states and
  tile opacity represents the within-species site share. Hybrid aggregation
  still uses rank transport and margin fitting, which endpoint plots do not
  reproduce.
