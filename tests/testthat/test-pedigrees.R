old_output <- options(msprimeR.output = "msprime_ts")

test_that("PedigreeBuilder builds pedigree tables", {
  pb <- PedigreeBuilder$new()
  mom <- pb$add_individual(time = 1)
  dad <- pb$add_individual(time = 1)
  kid <- pb$add_individual(time = 0, parents = c(mom, dad))
  expect_equal(c(mom, dad, kid), c(0, 1, 2))
  tables <- pb$finalise(sequence_length = 10)
  expect_equal(tables$sequence_length, 10)
  expect_equal(tables$individuals$parents, c(-1L, -1L, -1L, -1L, 0L, 1L))
  expect_equal(tables$nodes$flags, c(0, 0, 0, 0, 1, 1))
  expect_equal(tables$nodes$individual, rep(0:2, each = 2))
  expect_error(pb$add_individual(time = 0, parents = 1), "exactly two parents")
  multi <- PedigreeBuilder$new(demography = Demography$isolated_model(c(1, 1)))
  expect_error(multi$add_individual(time = 0), "Must specify the population")
})

test_that("pedigree file round trip and simulation", {
  f <- tempfile()
  writeLines(c("# id parent0 parent1 time is_sample",
               "a . . 2 0", "b . . 2 0", "c a b 1 0", "d a b 1 0",
               "e c d 0 1", "f c d 0 1"), f)
  tables <- parse_pedigree(f, sequence_length = 100)
  expect_equal(length(tables$individuals$flags), 6)
  off <- tables$individuals$metadata_offset
  md <- jsonlite::fromJSON(rawToChar(tables$individuals$metadata[(off[1] + 1):off[2]]))
  expect_equal(md$file_id, "a")
  out <- tempfile()
  write_pedigree(tables, out)
  again <- parse_pedigree(out, sequence_length = 100)
  expect_identical(again$individuals$parents, tables$individuals$parents)
  expect_identical(again$nodes, tables$nodes)
  ts <- sim_ancestry(initial_state = tables, model = "fixed_pedigree",
                     recombination_rate = 0.01, random_seed = 1)
  expect_equal(ts$num_samples, 4)

  bad <- tempfile()
  writeLines(c("# id parent0 parent1 time", "a x . 0"), bad)
  expect_error(parse_pedigree(bad), "Parent ID 'x' not defined")
  writeLines(c("# id parent0 parent1", "a . ."), bad)
  expect_error(parse_pedigree(bad), "columns are required")
  writeLines(c("# id parent0 parent1 time", "a . . 0", "a . . 0"), bad)
  expect_error(parse_pedigree(bad), "Duplicate ID")
})

test_that("fixed pedigree argument checks", {
  pb <- PedigreeBuilder$new()
  pb$add_individual(time = 0)
  tables <- pb$finalise(10)
  expect_error(sim_ancestry(initial_state = tables, model = "fixed_pedigree",
                            population_size = 10), "Cannot specify population_size")
  expect_error(sim_ancestry(2, model = "fixed_pedigree"), "input pedigree")
})

test_that("sim_pedigree", {
  fwd <- sim_pedigree(population_size = 10, end_time = 5, random_seed = 1,
                      sequence_length = 100)
  expect_equal(length(fwd$individuals$flags), 60)
  expect_equal(sum(fwd$nodes$flags == 1), 20)
  bwd <- sim_pedigree(population_size = 10, num_samples = 3, end_time = 4,
                      direction = "backward", random_seed = 1, sequence_length = 100)
  expect_equal(sum(bwd$nodes$flags == 1), 6)
  ts <- sim_ancestry(initial_state = fwd, model = "fixed_pedigree", random_seed = 2,
                     recombination_rate = 0.01)
  expect_equal(ts$num_samples, 20)
  set.seed(9)
  x <- stats::runif(1)
  set.seed(9)
  sim_pedigree(population_size = 5, end_time = 2, random_seed = 3)
  expect_equal(stats::runif(1), x)
  expect_error(sim_pedigree(population_size = 5, end_time = 2, direction = "sideways"),
               "unknown time direction")
})

options(old_output)
