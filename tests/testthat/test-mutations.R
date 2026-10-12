old_output <- options(msprimeR.output = "msprime_ts")

base_ts <- function() {
  sim_ancestry(5, population_size = 1e4, sequence_length = 1e4, recombination_rate = 1e-8,
               random_seed = 1)
}

test_that("matrix models are valid stochastic matrices", {
  models <- list(JC69(), HKY(2), F84(2), GTR(1:6), BLOSUM62(), PAM(), BinaryMutationModel(),
                 SMM(), TPM(p = 0.4, m = 0.5), EL2(m = 0.5, u = 0.4, v = 0.01),
                 MicrosatMutationModel(), InfiniteSites(NUCLEOTIDES),
                 JC69(state_independent = TRUE), F84(2, state_independent = TRUE),
                 GTR(1:6, state_independent = TRUE))
  for (m in models) {
    expect_equal(rowSums(m$transition_matrix), rep(1, length(m$alleles)), tolerance = 1e-12)
    expect_equal(sum(m$root_distribution), 1, tolerance = 1e-12)
    expect_true(all(m$transition_matrix >= 0))
  }
  expect_error(MatrixMutationModel(c("A", "B"), c(0.5, 0.6), diag(2)))
  expect_error(MatrixMutationModel(c("A", "B"), c(0.5, 0.5), diag(3)), "num_alleles")
  expect_error(TPM(p = 2, m = 0.5), "p must be between")
  expect_error(MicrosatMutationModel(lo = 5, hi = 2), "lo > hi")
  expect_error(InfiniteSites(3), "Bad alphabet")
})

test_that("sim_mutations basics", {
  ts <- base_ts()
  mts <- sim_mutations(ts, rate = 1e-7, random_seed = 2)
  expect_gt(mts$num_sites, 0)
  sites <- ts_table(mts, "sites")
  expect_true(all(sites$position == round(sites$position)))
  expect_true(all(sites$ancestral_state %in% c("A", "C", "G", "T")))
  cont <- sim_mutations(ts, rate = 1e-7, discrete_genome = FALSE, random_seed = 2)
  expect_false(all(ts_table(cont, "sites")$position == round(ts_table(cont, "sites")$position)))
  binary <- sim_mutations(ts, rate = 1e-7, model = "binary", random_seed = 2)
  expect_true(all(ts_table(binary, "sites")$ancestral_state == "0"))
  none <- sim_mutations(ts, random_seed = 2)
  expect_equal(none$num_sites, 0)
  late <- sim_mutations(ts, rate = 1e-6, start_time = 1e4, random_seed = 3)
  expect_true(all(ts_table(late, "mutations")$time >= 1e4))
  kept <- sim_mutations(mts, rate = 1e-7, random_seed = 4)
  expect_gt(kept$num_mutations, mts$num_mutations)
  replaced <- sim_mutations(mts, rate = 1e-7, keep = FALSE, random_seed = 4)
  expect_equal(replaced$num_provenances, 3)
})

test_that("sim_mutations argument checks", {
  ts <- base_ts()
  expect_error(sim_mutations(NULL), "TreeSequence")
  expect_error(sim_mutations(ts, rate = "a"), "float or a RateMap")
  expect_error(sim_mutations(ts, start_time = 10, end_time = 1), "start_time must be <=")
  expect_error(sim_mutations(ts, model = "nonsense"), "unknown")
  expect_error(sim_mutations(ts, rate = RateMap$uniform(10, 1)), "sequence length")
  expect_error(mutate(ts, model = JC69()), "legacy interface")
})

test_that("stateful SLiM models update next_id", {
  ts <- base_ts()
  m <- SLiMv6MutationModel()
  mts <- sim_mutations(ts, rate = 1e-7, model = m, random_seed = 3)
  expect_equal(m$next_id, mts$num_mutations)
  expect_warning(SLiMMutationModel(type = 1), "will be removed")
})

test_that("rate maps", {
  rm <- RateMap$new(position = c(0, 10, 20, 30), rate = c(0.1, NaN, 0.2))
  expect_equal(rm$sequence_length, 30)
  expect_equal(rm$num_intervals, 3)
  expect_equal(rm$num_missing_intervals, 1)
  expect_equal(rm$total_mass, 3)
  expect_equal(rm$mean_rate, 3 / 20)
  expect_equal(rm$get_rate(c(0, 15, 29.9)), c(0.1, NaN, 0.2))
  expect_equal(rm$get_cumulative_mass(c(0, 5, 10, 25, 30)), c(0, 0.5, 1, 2, 3))
  expect_equal(rm$find_index(15), 2)
  expect_equal(unname(rm$missing_intervals()), matrix(c(10, 20), 1))
  s <- rm$slice(5, 25)
  expect_equal(s$position, c(0, 5, 10, 20, 25, 30))
  expect_true(is.nan(s$rate[1]) && is.nan(s$rate[5]))
  st <- rm$slice(5, 25, trim = TRUE)
  expect_equal(st$position, c(0, 5, 15, 20))
  expect_error(RateMap$new(position = c(1, 2), rate = 1), "First position must be zero")
  expect_error(RateMap$new(position = c(0, 2, 1), rate = c(1, 1)), "strictly increasing")
  expect_error(RateMap$new(position = c(0, 1), rate = -1), "negative")
  expect_error(RateMap$new(position = c(0, 1), rate = NaN), "All intervals are missing")
  expect_output(print(rm), "RateMap")
})

test_that("read_hapmap", {
  f <- tempfile()
  writeLines(c("Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)",
               "chr1\t1000\t1\t0", "chr1\t2000\t2\t0.1", "chr1\t3000\t0\t0.3"), f)
  rm <- read_hapmap(f)
  expect_equal(rm$position, c(0, 1000, 2000, 3000))
  expect_true(is.nan(rm$rate[1]))
  expect_equal(rm$rate[2:3], c(1e-6, 2e-6))
  rm2 <- read_hapmap(f, rate_col = 2, sequence_length = 4000)
  expect_equal(rm2$sequence_length, 4000)
  expect_true(is.nan(rm2$rate[length(rm2$rate)]))
})

test_that("likelihoods behave", {
  arg <- sim_ancestry(3, population_size = 1, sequence_length = 5, recombination_rate = 0.1,
                      record_full_arg = TRUE, discrete_genome = FALSE, random_seed = 2)
  l1 <- log_arg_likelihood(arg, 0.1)
  l2 <- log_arg_likelihood(arg, 0.2)
  expect_true(is.finite(l1) && is.finite(l2))
  expect_error(log_arg_likelihood(arg, -1), "recombination_rate must be >= 0")
  mts <- sim_mutations(arg, rate = 0.5, discrete_genome = FALSE, random_seed = 3)
  expect_true(is.finite(log_mutation_likelihood(mts, 0.5)))
  if (mts$num_mutations > 0) expect_equal(log_mutation_likelihood(mts, 0), -Inf)
  expect_equal(log_mutation_likelihood(arg, 0), 0)
  plain <- sim_ancestry(3, population_size = 1, sequence_length = 5,
                        recombination_rate = 0.1, discrete_genome = FALSE, random_seed = 2)
  if (plain$num_trees > 1) expect_error(log_arg_likelihood(plain, 0.1), "NODE_IS_RE_EVENT")
})

options(old_output)
