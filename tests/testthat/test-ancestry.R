old_output <- options(msprimeR.output = "msprime_ts")

test_that("basic simulation produces a valid tree sequence", {
  ts <- sim_ancestry(10, population_size = 1e4, sequence_length = 1e5,
                     recombination_rate = 1e-8, random_seed = 1)
  expect_s3_class(ts, "msprime_ts")
  expect_equal(ts$num_samples, 20)
  expect_equal(ts$num_individuals, 10)
  expect_equal(ts$num_populations, 1)
  expect_equal(ts$sequence_length, 1e5)
  expect_equal(ts$time_units, "generations")
  expect_gt(ts$num_trees, 1)
  expect_equal(ts$num_provenances, 1)
  nodes <- ts_table(ts, "nodes")
  expect_equal(sum(nodes$flags == 1), 20)
  expect_true(all(nodes$time[nodes$flags == 1] == 0))
})

test_that("seeds give reproducible results, and set.seed controls unseeded runs", {
  a <- sim_ancestry(5, sequence_length = 100, recombination_rate = 0.01, random_seed = 5)
  b <- sim_ancestry(5, sequence_length = 100, recombination_rate = 0.01, random_seed = 5)
  expect_identical(a$tables$edges, b$tables$edges)
  set.seed(1)
  c1 <- sim_ancestry(5, sequence_length = 100, recombination_rate = 0.01)
  set.seed(1)
  c2 <- sim_ancestry(5, sequence_length = 100, recombination_rate = 0.01)
  expect_identical(c1$tables$edges, c2$tables$edges)
  expect_error(sim_ancestry(2, random_seed = 0), "seeds must be greater than 0")
  expect_error(sim_ancestry(2, random_seed = 2^32), "seeds must be")
})

test_that("sample specifications", {
  d <- Demography$island_model(c(100, 100), 0.1)
  ts1 <- sim_ancestry(c("0" = 2, "1" = 3), demography = d, random_seed = 2)
  pops <- ts_table(ts1, "nodes")$population[ts_table(ts1, "nodes")$flags == 1]
  expect_equal(as.vector(table(pops)), c(4, 6))
  d$populations[[1]]$name <- "A"
  d$populations[[2]]$name <- "B"
  ts2 <- sim_ancestry(c(A = 1, B = 2), demography = d, random_seed = 2)
  expect_equal(ts2$num_samples, 6)
  ts3 <- sim_ancestry(list(SampleSet(2, "A", time = 5), SampleSet(1, "B", ploidy = 3)),
                      demography = d, random_seed = 2)
  nodes <- ts_table(ts3, "nodes")
  samples <- nodes[nodes$flags == 1, ]
  expect_equal(nrow(samples), 7)
  expect_equal(sum(samples$time == 5), 4)
  expect_error(sim_ancestry(2, demography = d), "Numeric samples can only be used")
  expect_error(sim_ancestry(0, population_size = 1), "Zero samples specified")
  expect_error(sim_ancestry(list(SampleSet(2)), demography = d),
               "Must specify a SampleSet population")
  expect_error(sim_ancestry(c(Z = 1), demography = d), "not found")
})

test_that("argument validation mirrors msprime", {
  expect_error(sim_ancestry(), "Either the samples or initial_state")
  expect_error(sim_ancestry(2, recombination_rate = 1e-8), "sequence_length value must")
  expect_error(sim_ancestry(2, sequence_length = 10.5), "integer sequence length")
  expect_error(sim_ancestry(2, gene_conversion_rate = 1, sequence_length = 10),
               "tract length")
  expect_error(sim_ancestry(2, gene_conversion_tract_length = 1, sequence_length = 10),
               "gene conversion rate")
  expect_error(sim_ancestry(2, ploidy = 0), "ploidy must be >= 1")
  expect_error(sim_ancestry(2, ploidy = 1.5), "ploidy must be an integer")
  expect_error(sim_ancestry(2, model = "dtwf"), "population size must be set")
  expect_error(sim_ancestry(2, model = "nonsense"), "unknown")
  expect_error(sim_ancestry(2, record_full_arg = TRUE, coalescing_segments_only = FALSE),
               "Cannot set both")
  expect_error(sim_ancestry(2, additional_nodes = NodeType$RECOMBINANT),
               "coalescing_segments_only to FALSE")
  expect_error(sim_ancestry(2, demography = Demography$isolated_model(1),
                            population_size = 1), "Cannot specify demography and")
  expect_error(sim_ancestry(2, discrete_genome = "yes"), "Boolean flag")
  expect_error(sim_ancestry(2, model = SMCK(1), gene_conversion_rate = 1,
                            gene_conversion_tract_length = 1, sequence_length = 10),
               "not supported for the SMCK")
  expect_error(sim_ancestry(2, replicate_index = 1), "without a random_seed")
})

test_that("replicates", {
  reps <- sim_ancestry(3, num_replicates = 4, population_size = 10, random_seed = 3)
  expect_length(reps, 4)
  expect_true(all(vapply(reps, inherits, logical(1), "msprime_ts")))
  r2 <- sim_ancestry(3, replicate_index = 2, population_size = 10, random_seed = 3)
  expect_identical(r2$tables$edges, reps[[3]]$tables$edges)
  rec <- rawToChar(reps[[3]]$tables$provenances$record)
  expect_match(rec, "\"replicate_index\": 2")
  expect_length(sim_ancestry(3, num_replicates = 0, random_seed = 1), 0)
})

test_that("all ancestry models run", {
  models <- list("hudson", "smc", "smc_prime", SMCK(k = 10),
                 BetaCoalescent(alpha = 1.5), DiracCoalescent(psi = 0.5, c = 1),
                 DiscreteTimeWrightFisher())
  for (model in models) {
    ts <- sim_ancestry(4, population_size = 100, sequence_length = 100,
                       recombination_rate = 1e-3, model = model, random_seed = 7)
    expect_gte(ts$num_trees, 1)
  }
  ts <- sim_ancestry(4, population_size = 100, sequence_length = 100,
                     model = list(DiscreteTimeWrightFisher(duration = 10), "hudson"),
                     random_seed = 8)
  expect_equal(ts$num_samples, 8)
})

test_that("end_time and start_time", {
  ts <- sim_ancestry(5, population_size = 1000, end_time = 50, random_seed = 1)
  expect_lte(ts$max_time, 50)
  ts2 <- sim_ancestry(5, population_size = 1000, start_time = 100, random_seed = 1)
  nodes <- ts_table(ts2, "nodes")
  expect_true(all(nodes$time[nodes$flags != 1] > 100))
})

test_that("initial_state recapitation", {
  ts0 <- sim_ancestry(4, population_size = 100, end_time = 5, sequence_length = 100,
                      recombination_rate = 0.01, random_seed = 1)
  expect_gt(ts0$max_time, 0)
  ts <- sim_ancestry(initial_state = ts0, population_size = 100, random_seed = 2)
  expect_equal(ts$num_samples, ts0$num_samples)
  expect_error(sim_ancestry(initial_state = ts0, random_seed = 2),
               "Must specify either a demography")
  expect_error(sim_ancestry(2, initial_state = ts0, population_size = 1),
               "Cannot specify both samples and initial_state")
})

test_that("missing rate map flanks are trimmed", {
  rm <- RateMap$new(position = c(0, 10, 90, 100), rate = c(NaN, 0.01, NaN))
  ts <- sim_ancestry(3, population_size = 10, recombination_rate = rm, random_seed = 1)
  edges <- ts$tables$edges
  expect_true(all(edges$left >= 10 & edges$right <= 90))
  bad <- RateMap$new(position = c(0, 10, 20, 100), rate = c(0.1, NaN, 0.1))
  expect_error(sim_ancestry(3, population_size = 10, recombination_rate = bad),
               "Missing regions")
})

test_that("legacy simulate()", {
  ts <- simulate(6, Ne = 100, length = 1000, recombination_rate = 1e-4,
                 mutation_rate = 1e-4, random_seed = 1)
  expect_equal(ts$num_samples, 6)
  expect_gt(ts$num_sites, 0)
  expect_equal(ts$num_individuals, 0)
  reps <- simulate(4, num_replicates = 3, random_seed = 2)
  expect_length(reps, 3)
  expect_error(simulate(1), "Sample size must be >= 2")
  expect_error(simulate(), "Either sample_size")
  ts2 <- simulate(samples = list(Sample(0, 0), Sample(0, 1), c(0, 2)), random_seed = 3)
  expect_equal(ts2$num_samples, 3)
  ts3 <- simulate(4, demographic_events = list(SimulationModelChange(10, "smc")),
                  random_seed = 4)
  expect_equal(ts3$num_samples, 4)
  fit <- stats::lm(dist ~ speed, data = datasets::cars)
  expect_equal(dim(simulate(fit, nsim = 2, seed = 1)), c(50, 2))
  expect_error(simulate(4, nsim = 2), "unused argument")
})

test_that("census events record nodes", {
  d <- Demography$isolated_model(100)
  d$add_census(5)
  ts <- sim_ancestry(4, demography = d, sequence_length = 10, random_seed = 1,
                     coalescing_segments_only = FALSE, additional_nodes = NodeType$CENSUS)
  nodes <- ts_table(ts, "nodes")
  expect_true(any(nodes$time == 5 & nodes$flags == NODE_IS_CEN_EVENT))
})

test_that("full ARG recording flags nodes", {
  ts <- sim_ancestry(4, population_size = 10, sequence_length = 100,
                     recombination_rate = 0.01, record_full_arg = TRUE, random_seed = 3)
  flags <- ts$tables$nodes$flags
  expect_true(any(flags == NODE_IS_RE_EVENT))
})

test_that("migrations are recorded", {
  d <- Demography$island_model(c(10, 10), 0.5)
  ts <- sim_ancestry(c("0" = 2, "1" = 2), demography = d, record_migrations = TRUE,
                     random_seed = 1)
  expect_gt(ts$num_migrations, 0)
})

options(old_output)
