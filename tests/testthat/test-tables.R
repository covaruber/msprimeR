old_output <- options(msprimeR.output = "msprime_ts")

test_that("table dict round trips through C and .trees files", {
  ts <- sim_ancestry(4, population_size = 100, sequence_length = 100,
                     recombination_rate = 0.01, random_seed = 1)
  ts <- sim_mutations(ts, rate = 0.01, random_seed = 2)
  f <- tempfile(fileext = ".trees")
  ts_dump(ts, f)
  back <- ts_load(f)
  for (t in c("nodes", "edges", "sites", "mutations", "individuals", "populations",
              "provenances")) {
    expect_identical(back$tables[[t]], ts$tables[[t]])
  }
  expect_identical(ts_tables(f)$edges, ts$tables$edges)
  expect_error(ts_load(tempfile()), "reading tree sequence file")
})

test_that("ts_table decodes columns", {
  ts <- sim_mutations(sim_ancestry(3, sequence_length = 50, random_seed = 1),
                      rate = 0.1, random_seed = 1)
  sites <- ts_table(ts, "sites")
  expect_s3_class(sites, "data.frame")
  expect_equal(nrow(sites), ts$num_sites)
  expect_type(sites$ancestral_state, "character")
  prov <- ts_table(ts, "provenances")
  expect_equal(nrow(prov), 2)
  rec <- jsonlite::fromJSON(prov$record[1])
  expect_equal(rec$software$name, "msprimeR")
  expect_equal(rec$parameters$command, "sim_ancestry")
  expect_error(ts_table(ts, "nope"), "Unknown table")
})

test_that("invalid tables are rejected", {
  ts <- sim_ancestry(2, random_seed = 1)
  tables <- ts$tables
  tables$edges$parent[1] <- 1000L
  expect_error(msprimeR:::.new_ts(tables), "invalid tree sequence|out of bounds")
  tables <- ts$tables
  tables$nodes$time <- tables$nodes$time[-1]
  expect_error(msprimeR:::.new_ts(tables), "rows expected")
  expect_error(ts_tables(list(a = 1)), "Expected a tree sequence")
})

test_that("canonical JSON matches Python's encoding", {
  expect_equal(canonical_json(list(b = 1.5, a = "x")), "{\"a\":\"x\",\"b\":1.5}")
  expect_equal(canonical_json(list(z = TRUE, y = NULL, x = c(1, 2))),
               "{\"x\":[1,2],\"y\":null,\"z\":true}")
  expect_equal(canonical_json("caf\u00e9"), "\"caf\\u00e9\"")
  expect_equal(canonical_json(list(a = I(1))), "{\"a\":[1]}")
  expect_equal(canonical_json(1e-5), "1e-05")
})

test_that("tskitR interoperability", {
  skip_if_not_installed("tskitR")
  options(msprimeR.output = "tskitR")
  ts <- sim_ancestry(5, population_size = 100, sequence_length = 1000,
                     recombination_rate = 1e-4, random_seed = 1)
  expect_true(inherits(ts, "TreeSequence"))
  mts <- sim_mutations(ts, rate = 1e-4, random_seed = 2)
  expect_true(inherits(mts, "TreeSequence"))
  expect_gt(mts$num_sites, 0)
  m <- as_msprime_ts(mts)
  expect_true(as_tskit(m)$tables$equals(mts$tables))
  f <- tempfile(fileext = ".trees")
  mts$dump(f)
  expect_true(ts_load(f)$tables$equals(mts$tables))
  options(msprimeR.output = "msprime_ts")
})

test_that("tsdateR and tsinferR use the standard .trees handoff", {
  ts <- sim_mutations(
    sim_ancestry(12, population_size = 1000, sequence_length = 1000,
                 recombination_rate = 1e-6, random_seed = 4),
    rate = 1e-5, random_seed = 5)
  path <- tempfile(fileext = ".trees")
  ts_dump(ts, path)

  if (requireNamespace("tsdateR", quietly = TRUE)) {
    load_ts <- getExportedValue("tsdateR", "ts_load")
    date_ts <- getExportedValue("tsdateR", "date_ts")
    dump_ts <- getExportedValue("tsdateR", "ts_dump")
    dated <- date_ts(load_ts(path), mutation_rate = 1e-5, allow_unary = TRUE)
    expect_equal(getExportedValue("tsdateR", "ts_num_nodes")(dated), ts$num_nodes)
    dated_path <- tempfile(fileext = ".trees")
    dump_ts(dated, dated_path)
    if (requireNamespace("tskitR", quietly = TRUE)) {
      dated_tskit <- getExportedValue("tskitR", "tskit_load")(dated_path)
      expect_equal(dated_tskit$num_nodes, ts$num_nodes)
    }
  }

  if (requireNamespace("tsinferR", quietly = TRUE)) {
    load_ts <- getExportedValue("tsinferR", "ts_load")
    dump_ts <- getExportedValue("tsinferR", "ts_dump")
    info_ts <- getExportedValue("tsinferR", "ts_info")
    inferred_wrapper <- load_ts(path)
    expect_equal(info_ts(inferred_wrapper)$num_nodes, ts$num_nodes)
    returned_path <- tempfile(fileext = ".trees")
    dump_ts(inferred_wrapper, returned_path)
    expect_equal(msprimeR::ts_load(returned_path)$num_nodes, ts$num_nodes)
  }
})

test_that("output option is validated", {
  old <- options(msprimeR.output = "bad")
  on.exit(options(old))
  expect_error(sim_ancestry(2, random_seed = 1), "msprimeR.output")
})

test_that("library versions", {
  v <- library_versions()
  expect_named(v, c("gsl", "tskit", "kastore"))
})

options(old_output)
