# Build the parity fixtures used by tests/testthat/test-parity.R from the
# reference outputs of tools/parity/make_reference.py.
# Usage: Rscript tools/parity/make_fixtures.R <reference_dir> [lib]
args <- commandArgs(trailingOnly = TRUE)
ref <- args[1]
if (length(args) > 1) .libPaths(c(args[2], .libPaths()))
suppressPackageStartupMessages(library(msprimeR))
options(msprimeR.output = "msprime_ts")
source("tests/testthat/helper-scenarios.R")
names <- names(parity_scenarios())
digests <- vapply(names, function(nm) {
  parity_digest(file.path(ref, paste0(nm, ".trees")))
}, character(1))
dir.create("tests/testthat/fixtures", showWarnings = FALSE)
utils::write.csv(data.frame(scenario = names, digest = digests),
                 "tests/testthat/fixtures/parity_digests.csv", row.names = FALSE)
file.copy(file.path(ref, "values.json"), "tests/testthat/fixtures/parity_values.json",
          overwrite = TRUE)
cat("Wrote", length(digests), "digests\n")
