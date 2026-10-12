# msprimeR

msprimeR is an independent native R mirror of [msprime](https://tskit.dev/msprime/),
the coalescent and population genetics simulator. It bundles the original
msprime and tskit C libraries and does not depend on Python. This package is not
an official release of, or endorsed by, the upstream msprime project.

When reporting results produced with msprimeR, please cite the canonical
msprime publication, Baumdicker et al. (2022), *Efficient ancestry and mutation
simulation with msprime 1.0*, Genetics 220(3), iyab229,
<https://doi.org/10.1093/genetics/iyab229>. `citation("msprimeR")` provides
the full reference and the R mirror version.

* `sim_ancestry()`: coalescent with recombination and gene conversion,
  SMC/SMC'/SMC(k), Beta and Dirac multiple merger coalescents, selective sweeps,
  discrete-time Wright-Fisher and fixed pedigree simulations.
* `sim_mutations()`: JC69, HKY, F84, GTR, BLOSUM62, PAM, binary, microsatellite,
  infinite alleles and SLiM-compatible mutation models.
* `Demography`: populations, migration, splits, admixture, pulses, size
  changes; island and stepping stone models; species trees; demes models; the
  `DemographyDebugger`.
* `RateMap`, `PedigreeBuilder`, ARG and mutation likelihoods, and the legacy
  `simulate()`/`mutate()` interface.

Given the same parameters and random seed, msprimeR produces the same simulation
tables as upstream msprime. This package is an R mirror, not a separate
algorithmic publication; please cite the canonical msprime work above.

## Installation

msprimeR requires the GNU Scientific Library (GSL):

```sh
sudo apt-get install libgsl-dev   # Debian/Ubuntu
sudo dnf install gsl-devel        # Fedora
brew install gsl                  # macOS
```

then install the package from source, for example with
`R CMD INSTALL msprimeR` or `remotes::install_local("msprimeR")`.

## Example

```r
library(msprimeR)
ts <- sim_ancestry(samples = 10, population_size = 1e4, sequence_length = 1e6,
                   recombination_rate = 1e-8, random_seed = 1)
mts <- sim_mutations(ts, rate = 1e-8, random_seed = 2)
ts_dump(mts, "example.trees")
```

## Tree sequences and tskitR

Simulations return tree sequences. If the tskitR package is installed, these are
`tskitR::TreeSequence` objects (giving access to trees, statistics, genotypes,
and so on); otherwise, or with `options(msprimeR.output = "msprime_ts")`, a
lightweight `msprime_ts` object holding the tables is returned. Either can be
written to the standard `.trees` format with `ts_dump()`, which is readable by
tskitR, tsdateR, tsinferR, and other tskit-compatible R tools. See the
`interoperability` vignette for tested examples. tsinferR's inference pipeline
starts from genotype sources; it can read msprimeR outputs but does not infer
from an existing tree sequence.

## Differences from Python msprime

* Population, node and individual IDs are zero-based, as in the tree sequence
  tables. Interval indices returned by `RateMap$find_index()` are one-based.
* Classes with mutable state use R6: `Demography$new()`,
  `Demography$island_model(...)`, `demography$add_population(...)`,
  `RateMap$new(...)`, `PedigreeBuilder$new()`. Models are created with
  functions: `JC69()`, `SMCK(k = 10)`, `SampleSet(2, "A")`.
* Mappings are named vectors or lists: `samples = c(A = 2, B = 3)`.
* `num_replicates` returns a list rather than an iterator.
* Unseeded simulations draw their seed from R's RNG (`set.seed()` applies).

## Development

The C sources are vendored from the msprime and tskit repositories with
`tools/sync_c_sources.sh`. `tools/parity/` contains scripts that run the same
scenarios in Python msprime and msprimeR and compare the outputs table by table.
