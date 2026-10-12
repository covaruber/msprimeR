# Likelihood computations (port of msprime/likelihood.py).

#' Log-likelihood of mutations under infinite sites
#'
#' Returns the unnormalised log probability of the mutations on the tree
#' sequence under the infinite sites model with per-site, per-generation
#' mutation rate `mutation_rate`. Only the first mutation at each site is used.
#'
#' @param ts A tree sequence with mutations.
#' @param mutation_rate The non-negative mutation rate.
#' @return The log-likelihood (`-Inf` if the rate is zero and there are
#'   mutations).
#' @export
log_mutation_likelihood <- function(ts, mutation_rate) {
  .Call(msprimer_log_mutation_likelihood, ts_tables(ts), as.double(mutation_rate))
}

#' Log-likelihood of an ARG under the Hudson model
#'
#' Returns the log probability of the tree sequence under the coalescent with
#' recombination (Kuhner et al. 2000, eq. 1). The tree sequence must record the
#' full ARG (see `record_full_arg` in [sim_ancestry()]) on a continuous genome.
#'
#' @param ts A tree sequence recording the full ARG.
#' @param recombination_rate The per-link, per-generation recombination rate.
#' @param Ne The diploid effective population size.
#' @return The log-likelihood.
#' @export
log_arg_likelihood <- function(ts, recombination_rate, Ne = 1) {
  .Call(msprimer_log_likelihood_arg, ts_tables(ts), as.double(Ne),
        as.double(recombination_rate))
}
