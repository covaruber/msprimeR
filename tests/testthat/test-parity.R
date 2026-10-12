old_output <- options(msprimeR.output = "msprime_ts")

digests <- utils::read.csv(test_path("fixtures", "parity_digests.csv"),
                           stringsAsFactors = FALSE)
scenarios <- parity_scenarios()

test_that("all parity scenarios have reference digests", {
  expect_setequal(names(scenarios), digests$scenario)
})

for (k in seq_len(nrow(digests))) {
  local({
    name <- digests$scenario[k]
    expected <- digests$digest[k]
    test_that(paste("tables identical to Python msprime:", name), {
      skip_on_cran_if_slow <- name %in% c("dirac")
      if (skip_on_cran_if_slow) skip_on_cran()
      expect_identical(parity_digest(scenarios[[name]]()), expected)
    })
  })
}

values <- jsonlite::fromJSON(test_path("fixtures", "parity_values.json"),
                             simplifyVector = TRUE)

test_that("mutation model matrices match Python", {
  hky <- HKY(kappa = 2.0, equilibrium_frequencies = c(0.1, 0.2, 0.3, 0.4))
  expect_identical(unname(hky$transition_matrix), unname(values$hky_matrix))
  ms <- MicrosatMutationModel(s = 0.1, u = 0.5, v = 0.02, p = 0.8, m = 0.5, lo = 1, hi = 20)
  expect_equal(ms$root_distribution, values$microsat_root, tolerance = 1e-12)
})

test_that("stateful infinite alleles model continues numbering", {
  ia <- InfiniteAlleles()
  base <- parity_mutation_base()
  sim_mutations(base, rate = 1e-7, model = ia, random_seed = 32)
  sim_mutations(base, rate = 1e-7, model = ia, random_seed = 33)
  expect_equal(ia$next_allele, values$ia_next)
})

test_that("likelihoods match Python", {
  arg <- scenarios$lik_arg()
  expect_equal(log_arg_likelihood(arg, 0.1, Ne = 1), values$log_arg_likelihood,
               tolerance = 1e-12)
  mts <- sim_mutations(arg, rate = 0.5, discrete_genome = FALSE, random_seed = 41)
  expect_equal(log_mutation_likelihood(mts, 0.5), values$log_mutation_likelihood,
               tolerance = 1e-12)
})

test_that("demography debugger matches Python", {
  dbg <- parity_complex_demography()$debug()
  expect_equal(dbg$epoch_start_time, values$dbg_epoch_start)
  expect_equal(dbg$population_size_history, values$dbg_size_history, tolerance = 1e-12)
  expect_equal(dbg$population_size_trajectory(c(0, 10, 120, 500)), values$dbg_traj,
               tolerance = 1e-12)
  crt <- dbg$coalescence_rate_trajectory(seq(1, 1000, length.out = 30), c(A = 2, B = 2),
                                         double_step_validation = FALSE)
  expect_equal(crt$rates, values$dbg_rates, tolerance = 1e-9)
  expect_equal(crt$p, values$dbg_p, tolerance = 1e-9)
  P <- dbg$lineage_probabilities(c(10, 160, 350))
  expect_equal(P, values$dbg_lineage_probs, tolerance = 1e-9, ignore_attr = TRUE)
  dbg2 <- Demography$island_model(c(100, 200), 1e-2)$debug()
  expect_equal(dbg2$mean_coalescence_time(c("0" = 2)), values$mean_coal, tolerance = 1e-9)
})

options(old_output)
