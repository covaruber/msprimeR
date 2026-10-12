# Compare msprimeR outputs with reference outputs from Python msprime, table by
# table (see tools/parity/make_reference.py).
# Usage (from the package root): Rscript tools/parity/check_parity.R <reference_dir> [lib]
args <- commandArgs(trailingOnly = TRUE)
ref <- args[1]
if (length(args) > 1) .libPaths(c(args[2], .libPaths()))
suppressPackageStartupMessages(library(msprimeR))
options(msprimeR.output = "msprime_ts")
source("tests/testthat/helper-scenarios.R")

results <- list()
record <- function(name, ok, msg = "") {
  results[[name]] <<- ok
  cat(sprintf("%-22s %s %s\n", name, if (isTRUE(ok)) "OK  " else "FAIL", msg))
}

compare_tables <- function(name, ts) {
  py <- ts_tables(file.path(ref, paste0(name, ".trees")))
  r <- ts_tables(ts)
  diffs <- character(0)
  for (t in c("individuals", "nodes", "edges", "migrations", "sites", "mutations",
              "populations")) {
    if (!identical(py[[t]], r[[t]])) {
      cols <- names(py[[t]])[!mapply(identical, py[[t]], r[[t]][names(py[[t]])])]
      diffs <- c(diffs, sprintf("%s[%s]", t, paste(cols, collapse = ",")))
    }
  }
  for (f in c("sequence_length", "time_units")) {
    if (!identical(py[[f]], r[[f]])) diffs <- c(diffs, f)
  }
  record(name, length(diffs) == 0, paste(diffs, collapse = " "))
}

scenarios <- parity_scenarios()
for (nm in names(scenarios)) compare_tables(nm, scenarios[[nm]]())

reps <- sim_ancestry(3, population_size = 100, num_replicates = 3, sequence_length = 100,
                     recombination_rate = 1e-3, random_seed = 20)
for (j in 1:3) compare_tables(sprintf("rep_%d", j - 1), reps[[j]])

values <- jsonlite::fromJSON(file.path(ref, "values.json"), simplifyVector = TRUE)
ms <- MicrosatMutationModel(s = 0.1, u = 0.5, v = 0.02, p = 0.8, m = 0.5, lo = 1, hi = 20)
record("microsat_root", isTRUE(all.equal(ms$root_distribution, values$microsat_root,
                                         tolerance = 1e-12)))
hky <- HKY(kappa = 2.0, equilibrium_frequencies = c(0.1, 0.2, 0.3, 0.4))
record("hky_matrix", identical(unname(hky$transition_matrix), unname(values$hky_matrix)))
ia <- InfiniteAlleles()
base <- parity_mutation_base()
invisible(sim_mutations(base, rate = 1e-7, model = ia, random_seed = 32))
invisible(sim_mutations(base, rate = 1e-7, model = ia, random_seed = 33))
record("ia_next", ia$next_allele == values$ia_next)
arg <- scenarios$lik_arg()
record("log_arg_likelihood", isTRUE(all.equal(log_arg_likelihood(arg, 0.1, Ne = 1),
                                              values$log_arg_likelihood, tolerance = 1e-12)))
mts <- sim_mutations(arg, rate = 0.5, discrete_genome = FALSE, random_seed = 41)
compare_tables("lik_mut", mts)
record("log_mut_likelihood", isTRUE(all.equal(log_mutation_likelihood(mts, 0.5),
                                              values$log_mutation_likelihood,
                                              tolerance = 1e-12)))
dbg <- parity_complex_demography()$debug()
record("dbg_epoch_start", isTRUE(all.equal(dbg$epoch_start_time, values$dbg_epoch_start)))
record("dbg_size_history", isTRUE(all.equal(dbg$population_size_history,
                                            values$dbg_size_history, tolerance = 1e-12)))
record("dbg_traj", isTRUE(all.equal(dbg$population_size_trajectory(c(0, 10, 120, 500)),
                                    values$dbg_traj, tolerance = 1e-12)))
crt <- dbg$coalescence_rate_trajectory(seq(1, 1000, length.out = 30), c(A = 2, B = 2),
                                       double_step_validation = FALSE)
record("dbg_rates", isTRUE(all.equal(crt$rates, values$dbg_rates, tolerance = 1e-9)))
record("dbg_p", isTRUE(all.equal(crt$p, values$dbg_p, tolerance = 1e-9)))
P <- dbg$lineage_probabilities(c(10, 160, 350))
record("dbg_lineage_probs", isTRUE(all.equal(P, values$dbg_lineage_probs, tolerance = 1e-9,
                                             check.attributes = FALSE)))
dbg2 <- Demography$island_model(c(100, 200), 1e-2)$debug()
record("mean_coal", isTRUE(all.equal(dbg2$mean_coalescence_time(c("0" = 2)),
                                     values$mean_coal, tolerance = 1e-9)))

n_fail <- sum(!unlist(results))
cat(sprintf("\n%d checks, %d failures\n", length(results), n_fail))
quit(status = if (n_fail > 0) 1 else 0)
