# msprimeR 0.1.0

* Initial release: native R port of msprime 1.4.5 (development version),
  bundling the msprime C library and tskit C library 1.3.1.
* Issues found in msprime while porting, fixed here:
  - `log_mutation_likelihood()` could loop forever when a mutation lies below a
    unary root.
  - `DemographyDebugger$coalescence_rate_trajectory()` applied lineage movements
    at time zero twice when 0 was in `steps`, and used the coalescence rates of
    the previous step at the final time point.
  - `write_pedigree()` wrote `is_sample` as `True`/`False`, which
    `parse_pedigree()` cannot read; it is now written as 0/1.
  - `Demography$from_demes()` added redundant population parameter change
    events when only one of size or growth rate changed between epochs.
  - `Demography$add_population()` left the model in an invalid state after a
    failed call (e.g. a duplicate name).
  - The `sum_internal_gc_tract_lengths` simulator statistic is no longer
    truncated to an integer.
  - `mutate()` with an unsupported model now points to `sim_mutations()`.
  - The vendored `avl.c` no longer returns a value from a void function.
