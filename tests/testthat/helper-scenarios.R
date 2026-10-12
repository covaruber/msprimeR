# Scenarios shared by the parity tests (tests/testthat/test-parity.R) and the
# full comparison against Python msprime (tools/parity/check_parity.R). Each
# scenario mirrors one in tools/parity/make_reference.py and must produce
# tables identical to Python msprime.

parity_complex_demography <- function() {
  d <- Demography$new()
  d$add_population(name = "A", initial_size = 1000, growth_rate = 0.01)
  d$add_population(name = "B", initial_size = 2000)
  d$add_population(name = "C", initial_size = 500)
  d$add_population(name = "ADMIX", initial_size = 800)
  d$add_population(name = "ANC", initial_size = 5000)
  d$set_symmetric_migration_rate(c("A", "B"), 1e-3)
  d$set_migration_rate("C", "B", 5e-4)
  d$add_admixture(50, derived = "ADMIX", ancestral = c("A", "B"), proportions = c(0.3, 0.7))
  d$add_population_parameters_change(100, initial_size = 300, population = "A",
                                     growth_rate = 0)
  d$add_mass_migration(150, source = "C", dest = "B", proportion = 0.2)
  d$add_migration_rate_change(200, rate = 2e-3, source = "A", dest = "C")
  d$add_symmetric_migration_rate_change(250, c("A", "B"), 0)
  d$add_population_split(300, derived = c("A", "B", "C"), ancestral = "ANC")
  d$add_census(400)
  d
}

parity_pedigree <- function() {
  pb <- PedigreeBuilder$new()
  founders <- vapply(1:4, function(i) pb$add_individual(time = 3), numeric(1))
  g2 <- vapply(1:4, function(i) {
    pb$add_individual(time = 2, parents = c(founders[i], founders[(i %% 4) + 1]))
  }, numeric(1))
  g1 <- vapply(1:4, function(i) {
    pb$add_individual(time = 1, parents = c(g2[i], g2[((i + 1) %% 4) + 1]))
  }, numeric(1))
  for (i in 1:3) pb$add_individual(time = 0, parents = c(g1[i], g1[i + 1]))
  pb$finalise(sequence_length = 1000)
}

parity_mutation_base <- function() {
  sim_ancestry(6, population_size = 1e4, sequence_length = 1e4, recombination_rate = 1e-8,
               random_seed = 30)
}

parity_mutation_models <- function() {
  list(
    jc69 = JC69(), hky = HKY(kappa = 2.0, equilibrium_frequencies = c(0.1, 0.2, 0.3, 0.4)),
    f84 = F84(kappa = 3.0, equilibrium_frequencies = c(0.3, 0.2, 0.2, 0.3)),
    gtr = GTR(relative_rates = 1:6), blosum62 = BLOSUM62(), pam = PAM(),
    binary = BinaryMutationModel(), jc69_si = JC69(state_independent = TRUE),
    hky_si = HKY(kappa = 1.5, state_independent = TRUE), smm = SMM(lo = 2, hi = 10),
    tpm = TPM(p = 0.5, m = 0.3, lo = 2, hi = 10),
    el2 = EL2(m = 0.4, u = 0.6, v = 0.01, lo = 2, hi = 10),
    microsat = MicrosatMutationModel(s = 0.1, u = 0.5, v = 0.02, p = 0.8, m = 0.5,
                                     lo = 1, hi = 20),
    infinite_alleles = InfiniteAlleles(start_allele = 5),
    slim_v5 = SLiMv5MutationModel(type = 2, next_id = 10),
    slim_v6 = SLiMv6MutationModel(next_id = 3))
}

parity_demes_yaml <- "
time_units: years
generation_time: 25
demes:
  - name: ANC
    epochs: [{end_time: 5000, start_size: 2000}]
  - name: X
    ancestors: [ANC]
    epochs: [{start_size: 1000, end_size: 4000}]
  - name: Y
    ancestors: [ANC]
    epochs: [{end_time: 1000, start_size: 500}, {start_size: 800}]
migrations:
  - demes: [X, Y]
    rate: 1e-4
pulses:
  - {sources: [X], dest: Y, time: 2000, proportions: [0.1]}
"

# Named list of functions, each returning a tree sequence.
parity_scenarios <- function() {
  s <- list(
    basic = function() sim_ancestry(5, sequence_length = 1e5, recombination_rate = 1e-8,
                                    population_size = 1e4, random_seed = 42),
    island = function() {
      sim_ancestry(c("0" = 3, "1" = 2, "2" = 2),
                   demography = Demography$island_model(c(1000, 1000, 500), 1e-3),
                   sequence_length = 5e4, recombination_rate = 1e-8, record_migrations = TRUE,
                   discrete_genome = FALSE, random_seed = 7)
    },
    complex = function() {
      sim_ancestry(list(SampleSet(3, "A"), SampleSet(2, "B", time = 10),
                        SampleSet(2, "C", ploidy = 1), SampleSet(2, "ADMIX")),
                   demography = parity_complex_demography(), sequence_length = 1e5,
                   recombination_rate = 1e-8, random_seed = 11)
    },
    dtwf = function() sim_ancestry(4, population_size = 100, model = "dtwf",
                                   sequence_length = 1e4, recombination_rate = 1e-6,
                                   random_seed = 3),
    multi_model = function() {
      sim_ancestry(4, population_size = 100, sequence_length = 1e4, recombination_rate = 1e-6,
                   model = list(DiscreteTimeWrightFisher(duration = 20), StandardCoalescent()),
                   random_seed = 4)
    },
    smc = function() sim_ancestry(5, population_size = 1e4, model = "smc",
                                  sequence_length = 1e5, recombination_rate = 1e-8,
                                  random_seed = 5),
    smc_prime = function() sim_ancestry(5, population_size = 1e4, model = "smc_prime",
                                        sequence_length = 1e5, recombination_rate = 1e-8,
                                        random_seed = 5),
    smck = function() sim_ancestry(5, population_size = 1e4, model = SMCK(k = 100),
                                   sequence_length = 1e5, recombination_rate = 1e-8,
                                   random_seed = 6),
    beta = function() sim_ancestry(6, population_size = 1e4,
                                   model = BetaCoalescent(alpha = 1.5),
                                   sequence_length = 1e4, recombination_rate = 1e-8,
                                   random_seed = 8),
    dirac = function() sim_ancestry(6, population_size = 1e4,
                                    model = DiracCoalescent(psi = 0.5, c = 1),
                                    sequence_length = 1e4, recombination_rate = 1e-8,
                                    random_seed = 9),
    sweep = function() {
      sim_ancestry(5, population_size = 1e3, sequence_length = 1e4, recombination_rate = 1e-7,
                   model = list(SweepGenicSelection(position = 5000, start_frequency = 1 / 2000,
                                                    end_frequency = 0.99, s = 0.25, dt = 1e-6),
                                StandardCoalescent()), random_seed = 10)
    },
    gc = function() sim_ancestry(5, population_size = 1e4, sequence_length = 1e5,
                                 recombination_rate = 1e-8, gene_conversion_rate = 1e-8,
                                 gene_conversion_tract_length = 100, random_seed = 12),
    full_arg = function() sim_ancestry(4, population_size = 1e4, sequence_length = 1e4,
                                       recombination_rate = 1e-7, record_full_arg = TRUE,
                                       discrete_genome = FALSE, random_seed = 13),
    additional_nodes = function() {
      sim_ancestry(4, population_size = 1e4, sequence_length = 1e4, recombination_rate = 1e-7,
                   coalescing_segments_only = FALSE,
                   additional_nodes = NodeType$RECOMBINANT + NodeType$COMMON_ANCESTOR,
                   random_seed = 14)
    },
    end_time = function() sim_ancestry(6, population_size = 1e4, end_time = 5000,
                                       sequence_length = 1e4, random_seed = 15),
    ploidy1 = function() sim_ancestry(6, ploidy = 1, population_size = 100, random_seed = 16),
    ploidy4 = function() sim_ancestry(3, ploidy = 4, population_size = 100,
                                      sequence_length = 1000, recombination_rate = 1e-4,
                                      random_seed = 17),
    no_local_mrca = function() sim_ancestry(4, population_size = 100, sequence_length = 1000,
                                            recombination_rate = 1e-4, random_seed = 18,
                                            stop_at_local_mrca = FALSE),
    ratemap = function() {
      rm <- RateMap$new(position = c(0, 1000, 5000, 9000, 10000),
                        rate = c(NaN, 1e-7, 5e-8, NaN))
      sim_ancestry(4, population_size = 1e4, recombination_rate = rm, random_seed = 19)
    },
    rep_index = function() sim_ancestry(3, population_size = 100, replicate_index = 2,
                                        sequence_length = 100, recombination_rate = 1e-3,
                                        random_seed = 20),
    recap = function() {
      ts0 <- sim_ancestry(4, population_size = 100, end_time = 10, sequence_length = 1000,
                          recombination_rate = 1e-4, random_seed = 21)
      sim_ancestry(initial_state = ts0, population_size = 200, recombination_rate = 1e-4,
                   random_seed = 22)
    },
    pedigree_tables = function() as_msprime_ts(parity_pedigree()),
    pedigree_sim = function() sim_ancestry(initial_state = parity_pedigree(),
                                           model = "fixed_pedigree",
                                           recombination_rate = 1e-3, random_seed = 23),
    simulate_basic = function() simulate(10, Ne = 1000, length = 1e4,
                                         recombination_rate = 1e-7, mutation_rate = 1e-7,
                                         random_seed = 24),
    simulate_popconf = function() {
      simulate(population_configurations = list(
                 PopulationConfiguration(sample_size = 4, initial_size = 500),
                 PopulationConfiguration(sample_size = 3, initial_size = 500,
                                         growth_rate = 0.01)),
               migration_matrix = matrix(c(0, 1e-3, 1e-3, 0), 2),
               demographic_events = list(
                 MassMigration(time = 200, source = 1, dest = 0, proportion = 1),
                 CensusEvent(time = 300)),
               length = 1e4, recombination_rate = 1e-7, random_seed = 25)
    },
    mut_base = parity_mutation_base,
    mut_cont = function() sim_mutations(parity_mutation_base(), rate = 1e-7,
                                        discrete_genome = FALSE, random_seed = 34),
    mut_ratemap = function() {
      mrm <- RateMap$new(position = c(0, 3000, 10000), rate = c(1e-7, 5e-7))
      sim_mutations(parity_mutation_base(), rate = mrm, start_time = 100, end_time = 20000,
                    random_seed = 35)
    },
    mut_keep = function() {
      sim_mutations(sim_mutations(parity_mutation_base(), rate = 1e-7, random_seed = 36),
                    rate = 1e-7, random_seed = 37)
    },
    mutate_legacy = function() mutate(parity_mutation_base(), rate = 1e-7, random_seed = 38),
    mut_ia_chain = function() {
      ia <- InfiniteAlleles()
      base <- parity_mutation_base()
      sim_mutations(base, rate = 1e-7, model = ia, random_seed = 32)
      sim_mutations(base, rate = 1e-7, model = ia, random_seed = 33)
    },
    lik_arg = function() sim_ancestry(4, population_size = 1, sequence_length = 10,
                                      recombination_rate = 0.1, record_full_arg = TRUE,
                                      discrete_genome = FALSE, random_seed = 40),
    species_tree = function() {
      st <- Demography$from_species_tree("((A:10.0,B:10.0)AB:20.0,C:30.0)", initial_size = 100)
      sim_ancestry(c(A = 2, B = 2, C = 2), demography = st, random_seed = 50)
    },
    demes = function() {
      dm <- Demography$from_demes(parity_demes_yaml)
      sim_ancestry(c(X = 2, Y = 2), demography = dm, sequence_length = 1e4,
                   recombination_rate = 1e-8, random_seed = 51)
    }
  )
  for (nm in names(parity_mutation_models())) {
    local({
      model_name <- nm
      s[[paste0("mut_", model_name)]] <<- function() {
        sim_mutations(parity_mutation_base(), rate = 1e-7,
                      model = parity_mutation_models()[[model_name]], random_seed = 31)
      }
    })
  }
  s
}

# Digest of the simulation-relevant tables (provenance excluded), computed from
# the raw binary representation of each column.
parity_digest <- function(ts) {
  tables <- ts_tables(ts)
  f <- tempfile()
  on.exit(unlink(f))
  con <- file(f, "wb")
  for (t in c("individuals", "nodes", "edges", "migrations", "sites", "mutations",
              "populations")) {
    tab <- tables[[t]]
    for (col in sort(names(tab))) {
      v <- tab[[col]]
      if (is.character(v)) v <- charToRaw(enc2utf8(paste(v, collapse = "\n")))
      writeBin(c(t, col), con)
      writeBin(v, con, endian = "little")
    }
  }
  writeBin(c(as.character(tables$sequence_length), tables$time_units), con)
  close(con)
  unname(tools::md5sum(f))
}
