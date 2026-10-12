test_that("Demography construction and lookup", {
  d <- Demography$new()
  a <- d$add_population(name = "A", initial_size = 100, description = "first")
  d$add_population(initial_size = 200)
  expect_equal(d$num_populations, 2)
  expect_equal(names(d), c("A", "pop_1"))
  expect_identical(d$get("A"), a)
  expect_identical(d[["A"]], a)
  expect_identical(d[[1]], d$get("pop_1"))
  expect_equal(d$get(0)$description, "first")
  d[["A"]]$initial_size <- 50
  expect_equal(d$populations[[1]]$initial_size, 50)
  expect_error(d$get("Z"), "not found")
  expect_error(d$get(2), "out of bounds")
  expect_error(d$add_population(name = "A", initial_size = 1), "Duplicate population name")
  expect_error(d$add_population(name = "1bad", initial_size = 1), "valid identifier")
  expect_error(d$add_population(initial_size = -1), "Negative population size")
  # Failed additions leave the model unchanged.
  expect_equal(d$num_populations, 2)
  expect_equal(dim(d$migration_matrix), c(2, 2))
})

test_that("R6 basics still work with the [[ method", {
  d <- Demography$isolated_model(c(1, 2))
  expect_output(print(d), "Demography")
  expect_true(is.function(d[["add_population"]]))
  expect_length(d, 2)
})

test_that("migration rates", {
  d <- Demography$isolated_model(c(1, 1, 1))
  d$set_migration_rate(0, 1, 0.5)
  expect_equal(d$migration_matrix[1, 2], 0.5)
  d$set_symmetric_migration_rate(c(0, 1, 2), 0.1)
  expect_equal(d$migration_matrix, matrix(c(0, .1, .1, .1, 0, .1, .1, .1, 0), 3))
  expect_error(d$set_migration_rate(0, 0, 1), "must be different")
  expect_error(d$set_symmetric_migration_rate(0, 1), "at least two")
  island <- Demography$island_model(c(1, 1, 1), 0.2)
  expect_equal(sum(island$migration_matrix), 6 * 0.2)
  ss <- Demography$stepping_stone_model(rep(1, 4), 0.1)
  expect_equal(sum(ss$migration_matrix > 0), 8)
  ssb <- Demography$stepping_stone_model(rep(1, 4), 0.1, boundaries = TRUE)
  expect_equal(sum(ssb$migration_matrix > 0), 6)
})

test_that("events and validation", {
  d <- Demography$isolated_model(c(1, 1, 1))
  d$add_population_split(10, derived = c(0, 1), ancestral = 2)
  expect_false(d$get(2)$initially_active)
  expect_equal(d$get(2)$default_sampling_time, 10)
  d$add_mass_migration(5, source = 0, dest = 1, proportion = 0.5)
  expect_error(d$validate(), "time-sorted")
  d$sort_events()
  r <- d$validate()
  expect_equal(r$get(0)$default_sampling_time, 0)
  expect_true(r$get(0)$initially_active)
  expect_null(d$get(0)$default_sampling_time)
  expect_error(d$add_admixture(1, derived = 0, ancestral = c(1, 2), proportions = c(0.5, 0.4)),
               "approximately 1")
  expect_error(d$add_population_parameters_change(1), "one or more")
  expect_error(d$add_mass_migration(1, source = 0, dest = 7, proportion = 1), "out of bounds")
  expect_output(print(d), "Population Split")
})

test_that("copy reorders populations and remaps events", {
  d <- Demography$new()
  d$add_population(name = "A", initial_size = 1)
  d$add_population(name = "B", initial_size = 2)
  d$add_population(name = "C", initial_size = 3)
  d$set_migration_rate("A", "B", 0.1)
  d$add_population_split(10, derived = c("A", "B"), ancestral = "C")
  c2 <- d$copy(c("C", "B", "A"))
  expect_equal(names(c2), c("C", "B", "A"))
  expect_equal(c2$migration_matrix[3, 2], 0.1)
  expect_true(d$is_equivalent(d$copy()))
  expect_true(d$equals(d$copy()))
})

test_that("equivalence checking", {
  d1 <- Demography$isolated_model(c(10, 10))
  d1$add_population_parameters_change(5, initial_size = 3, population = 0)
  d2 <- Demography$isolated_model(c(10, 10))
  d2$add_population_parameters_change(5, initial_size = 3, population = 0)
  expect_true(d1$is_equivalent(d2))
  d2$add_population_parameters_change(6, initial_size = 4, population = 0)
  expect_false(d1$is_equivalent(d2))
})

test_that("debugger epochs", {
  d <- Demography$isolated_model(c(100, 200), growth_rate = c(0.01, 0))
  d$add_population_parameters_change(10, initial_size = 50, growth_rate = 0, population = 0)
  dbg <- d$debug()
  expect_equal(dbg$num_epochs, 2)
  expect_equal(dbg$epoch_start_time, c(0, 10))
  expect_equal(dbg$epochs[[1]]$populations[[1]]$end_size, 100 * exp(-0.01 * 10))
  expect_equal(dbg$population_size_trajectory(c(0, 5))[, 1], 100 * exp(-0.01 * c(0, 5)))
  expect_output(print(dbg), "Epoch\\[1\\]")
  locs <- dbg$possible_lineage_locations()
  expect_true(all(locs$indicators))
})

test_that("misspecification warning after merge", {
  d <- Demography$isolated_model(c(1, 1))
  d$set_migration_rate(0, 1, 0.1)
  d$add_mass_migration(5, source = 0, dest = 1, proportion = 1)
  expect_warning(d$debug(), "Non-zero migration rates exist after merging")
})

test_that("from_tree_sequence recovers population names", {
  d <- Demography$new()
  d$add_population(name = "north", initial_size = 10, description = "n")
  d$add_population(name = "south", initial_size = 10)
  d$set_symmetric_migration_rate(c("north", "south"), 0.1)
  ts <- sim_ancestry(c(north = 1, south = 1), demography = d, random_seed = 1,
                     record_provenance = FALSE)
  d2 <- Demography$from_tree_sequence(ts)
  expect_equal(names(d2), c("north", "south"))
  expect_equal(d2$get("north")$description, "n")
})

test_that("extra metadata is stored in the population table", {
  d <- Demography$new()
  d$add_population(name = "A", initial_size = 10, extra_metadata = list(x = 1, y = "z"))
  ts <- sim_ancestry(2, demography = d, random_seed = 1)
  md <- rawToChar(ts_tables(ts)$populations$metadata)
  expect_equal(md, "{\"description\":\"\",\"name\":\"A\",\"x\":1,\"y\":\"z\"}")
  d$add_population(name = "B", initial_size = 1, extra_metadata = list(name = "B"))
  expect_error(sim_ancestry(c(A = 1), demography = d), "Cannot set standard metadata")
})

test_that("from_old_style", {
  pcs <- list(PopulationConfiguration(initial_size = 10, metadata = list(name = "X")),
              PopulationConfiguration(initial_size = 20))
  d <- Demography$from_old_style(pcs, migration_matrix = matrix(c(0, 1, 1, 0), 2))
  expect_equal(names(d), c("X", "pop_1"))
  expect_error(Demography$from_old_style(list(PopulationConfiguration(sample_size = 2))),
               "sample_size")
})

test_that("species trees", {
  d <- Demography$from_species_tree("(A:10,B:10)C;", initial_size = 100)
  expect_equal(names(d), c("A", "B", "C"))
  expect_equal(d$events[[1]]$time, 10)
  d2 <- Demography$from_species_tree("(A:1,B:1)", initial_size = c(A = 1, B = 2, .default = 3),
                                     time_units = "myr", generation_time = 25)
  expect_equal(d2$events[[1]]$time, 1e6 / 25)
  expect_equal(d2$get(2)$initial_size, 3)
  expect_error(Demography$from_species_tree("(A:1,B:2)", 1), "time 0")
  expect_error(Demography$from_species_tree("(A:1,B:1)", 1, time_units = "x"), "not accepted")
  expect_error(Demography$from_species_tree("(A:1,B:1)", 1, time_units = "yr"),
               "generation time must be")
})

test_that("StarBEAST nexus trees", {
  nexus <- paste(
    "#NEXUS", "Begin trees;", "Translate", "1 spc1,", "2 spc2,", "3 spc3", ";",
    "tree TREE1 = ((1[&dmv={0.1}]:1.0,2[&dmv={0.2}]:1.0)[&dmv={0.3}]:1.0,3[&dmv={0.4}]:2.0)[&dmv={0.5}];",
    "End;", sep = "\n")
  d <- Demography$from_starbeast(nexus, generation_time = 5)
  expect_equal(names(d)[1:3], c("spc1", "spc2", "spc3"))
  expect_equal(d$get("spc1")$initial_size, 0.1 * 1e6 / 5)
  expect_equal(d$num_events, 2)
})

test_that("demes round trip", {
  skip_if_not_installed("yaml")
  d <- Demography$new()
  d$add_population(name = "A", initial_size = 100)
  d$add_population(name = "B", initial_size = 200)
  d$add_population(name = "C", initial_size = 300)
  d$set_symmetric_migration_rate(c("A", "B"), 1e-3)
  d$add_population_split(100, derived = c("A", "B"), ancestral = "C")
  g <- d$to_demes()
  expect_s3_class(g, "demes_graph")
  expect_length(g$demes, 3)
  d2 <- Demography$from_demes(g)
  # demes orders demes so that ancestors come first.
  expect_equal(names(d2), c("C", "A", "B"))
  expect_true(d$is_equivalent(d2$copy(names(d))))
  f <- tempfile(fileext = ".yaml")
  write_demes(g, f)
  d3 <- Demography$from_demes(f)
  expect_true(d$is_equivalent(d3$copy(names(d))))
  expect_error(Demography$from_demes(list(time_units = "years", demes = list())),
               "generation_time must be specified")
})

test_that("demes yaml strings with YAML 1.1 pitfalls", {
  skip_if_not_installed("yaml")
  y <- "time_units: generations\ndemes:\n  - name: Y\n    epochs: [{start_size: 1e3}]\n"
  d <- Demography$from_demes(y)
  expect_equal(names(d), "Y")
  expect_equal(d$get("Y")$initial_size, 1000)
})
