# The deprecated msprime 0.x simulate() interface (port of
# msprime.ancestry.simulate) and the legacy RecombinationMap.

#' Recombination map (deprecated)
#'
#' Legacy representation of a recombination map, for use with [simulate()].
#' `rates` has the same length as `positions` and its last value is ignored.
#' Please use [RateMap] in new code.
#'
#' @param positions The positions where the rate changes (starting at 0).
#' @param rates The rates (same length as positions; the last is ignored).
#' @param num_loci Must be NULL or equal to the last position (in which case
#'   a discrete genome is simulated).
#' @return A `RecombinationMap` object, a list holding a [RateMap] in `map`.
#' @export
RecombinationMap <- function(positions, rates, num_loci = NULL) {
  is_discrete <- !is.null(num_loci) && num_loci == positions[length(positions)]
  if (!is.null(num_loci) && num_loci != positions[length(positions)]) {
    stop("The RecombinationMap interface is deprecated and only partially supported. If ",
         "you wish to simulate a number of discrete loci, you must set num_loci == the ",
         "sequence length. Otherwise, num_loci is no longer supported.", call. = FALSE)
  }
  structure(list(map = RateMap$new(position = positions, rate = rates[-length(rates)]),
                 is_discrete = is_discrete),
            class = c("RecombinationMap", "msprimeR_object"))
}

#' @rdname RecombinationMap
#' @param length The sequence length.
#' @param rate The uniform recombination rate.
#' @export
uniform_recombination_map <- function(length, rate, num_loci = NULL) {
  RecombinationMap(c(0, length), c(rate, 0), num_loci = num_loci)
}

.samples_factory <- function(sample_size, samples, population_configurations) {
  the_samples <- list()
  if (!is.null(sample_size)) {
    if (!is.null(samples)) {
      stop("Cannot specify sample size and samples simultaneously.", call. = FALSE)
    }
    if (!is.null(population_configurations)) {
      stop("Cannot specify sample size and population_configurations simultaneously.",
           call. = FALSE)
    }
    the_samples <- rep(list(.Sample(0, 0)), sample_size)
  }
  if (!is.null(population_configurations)) {
    for (pc in population_configurations) {
      if (!inherits(pc, "PopulationConfiguration")) {
        stop("Population configurations must be a list of PopulationConfiguration ",
             "instances", call. = FALSE)
      }
    }
    if (is.null(samples)) {
      the_samples <- list()
      for (j in seq_along(population_configurations)) {
        ss <- population_configurations[[j]]$sample_size
        if (!is.null(ss)) the_samples <- c(the_samples, rep(list(.Sample(j - 1, 0)), ss))
      }
    } else {
      for (pc in population_configurations) {
        if (!is.null(pc$sample_size)) {
          stop("Cannot specify population configuration sample size and samples ",
               "simultaneously", call. = FALSE)
        }
      }
      the_samples <- samples
    }
  } else if (!is.null(samples)) {
    the_samples <- samples
  }
  lapply(the_samples, function(s) {
    if (inherits(s, "Sample")) return(s)
    if (length(s) == 2) return(.Sample(s[[1]], s[[2]]))
    stop("Each sample must be a Sample or a (population, time) pair", call. = FALSE)
  })
}

.build_initial_tables_simulate <- function(sequence_length, samples, demography) {
  tables <- .empty_tables(sequence_length)
  tables$time_units <- TIME_UNITS_GENERATIONS
  for (index in seq_along(samples)) {
    pop <- samples[[index]]$population
    if (pop < 0) .stopf("Negative population ID in sample at index %d", index - 1)
    if (pop >= demography$num_populations) {
      .stopf("Invalid population reference '%s' in sample at index %d", pop, index - 1)
    }
  }
  n <- length(samples)
  tables$nodes <- .append_rows(tables$nodes, list(
    flags = rep(1, n),
    time = vapply(samples, function(s) as.double(s$time), numeric(1)),
    population = vapply(samples, function(s) as.integer(s$population), integer(1)),
    individual = rep(-1L, n)))
  md <- lapply(demography$populations, function(p) {
    if (is.null(p$extra_metadata)) raw(0) else
      charToRaw(.json_encode(.json_object(p$extra_metadata), sort_keys = FALSE,
                             item_sep = ", ", key_sep = ": "))
  })
  tables$populations <- .append_rows(tables$populations, list(), list(metadata = md))
  tables
}

#' Simulate the coalescent (deprecated 0.x interface)
#'
#' The msprime 0.x `simulate()` function, supported for compatibility. Note
#' that `Ne` is the diploid effective population size and `sample_size` the
#' number of monoploid genomes. Please use [sim_ancestry()] and
#' [sim_mutations()] in new code.
#'
#' @param sample_size The number of sampled monoploid genomes.
#' @param Ne The effective (diploid) population size (default 1).
#' @param length The sequence length (default 1).
#' @param recombination_rate The recombination rate per unit length.
#' @param recombination_map A [RateMap] or [RecombinationMap()].
#' @param mutation_rate If specified, binary infinite sites mutations are
#'   added at this rate.
#' @param population_configurations A list of [PopulationConfiguration()].
#' @param migration_matrix The migration matrix.
#' @param demographic_events A list of [demographic_events] (and
#'   [SimulationModelChange()] objects).
#' @param samples A list of [Sample()] objects (or (population, time) pairs).
#' @param model The ancestry model.
#' @param record_migrations If TRUE, record migrations.
#' @param random_seed The random seed.
#' @param replicate_index Return only this replicate.
#' @param mutation_generator No longer supported (must be NULL).
#' @param num_replicates If given, return a list of replicates.
#' @param from_ts A tree sequence to start from.
#' @param start_time,end_time Start and end times.
#' @param record_full_arg If TRUE, record the full ARG.
#' @param num_labels Number of labels (internal).
#' @param record_provenance If TRUE, record provenance.
#' @param ... If `sample_size` is a fitted model object (anything with a
#'   class), the call is passed on to [stats::simulate()] with these
#'   arguments, so that attaching msprimeR does not mask it. Otherwise unused.
#' @return A tree sequence, or a list of them.
#' @export
simulate <- function(sample_size = NULL, Ne = 1, length = NULL, recombination_rate = NULL,
                     recombination_map = NULL, mutation_rate = NULL,
                     population_configurations = NULL, migration_matrix = NULL,
                     demographic_events = NULL, samples = NULL, model = NULL,
                     record_migrations = FALSE, random_seed = NULL, replicate_index = NULL,
                     mutation_generator = NULL, num_replicates = NULL, from_ts = NULL,
                     start_time = NULL, end_time = NULL, record_full_arg = FALSE,
                     num_labels = NULL, record_provenance = TRUE, ...) {
  if (is.object(sample_size)) return(stats::simulate(sample_size, ...))
  if (...length() > 0) {
    stop("unused argument(s): ", paste(names(list(...)), collapse = ", "), call. = FALSE)
  }
  replicate_index <- .parse_replicate_index(replicate_index, random_seed, num_replicates)
  random_seed <- .parse_random_seed(random_seed)
  provenance_dict <- NULL
  if (record_provenance) {
    provenance_dict <- .get_provenance_dict(.parameters_dict("simulate", list(
      sample_size = sample_size, Ne = Ne, length = length,
      recombination_rate = recombination_rate, recombination_map = recombination_map,
      mutation_rate = mutation_rate, population_configurations = population_configurations,
      migration_matrix = migration_matrix, demographic_events = demographic_events,
      samples = samples, model = model, record_migrations = record_migrations,
      from_ts = from_ts, start_time = start_time, end_time = end_time,
      record_full_arg = record_full_arg, num_labels = num_labels,
      random_seed = random_seed)))
  }
  if (!is.null(mutation_generator)) {
    stop("mutation_generator is not longer supported. Please use mutate instead",
         call. = FALSE)
  }
  if (!is.null(mutation_rate)) {
    if (!is.null(from_ts)) {
      stop("Cannot specify mutation rate combined with from_ts. Please use mutate on the ",
           "final tree sequence instead", call. = FALSE)
    }
    if (!is.null(start_time) && start_time > 0) {
      stop("Cannot specify mutation rate combined with a non-zero start_time. Please use ",
           "mutate on the returned tree sequence instead", call. = FALSE)
    }
    mutation_rate <- as.double(mutation_rate)
  }
  sim <- .parse_simulate(
    sample_size = sample_size, Ne = Ne, length = length,
    recombination_rate = recombination_rate, recombination_map = recombination_map,
    population_configurations = population_configurations,
    migration_matrix = migration_matrix, demographic_events = demographic_events,
    samples = samples, model = model, record_migrations = record_migrations,
    from_ts = from_ts, start_time = start_time, end_time = end_time,
    record_full_arg = record_full_arg, num_labels = num_labels, random_seed = random_seed)
  on.exit(sim$free(), add = TRUE)
  .wrap_replicates(sim, num_replicates, replicate_index, provenance_dict,
                   mutation_rate = mutation_rate)
}

.parse_simulate <- function(sample_size = NULL, Ne = 1, length = NULL,
                            recombination_rate = NULL, recombination_map = NULL,
                            population_configurations = NULL, migration_matrix = NULL,
                            samples = NULL, demographic_events = NULL, model = NULL,
                            record_migrations = FALSE, from_ts = NULL, start_time = NULL,
                            end_time = NULL, record_full_arg = FALSE, num_labels = NULL,
                            random_seed = NULL) {
  if (Ne <= 0) stop("Population size must be positive", call. = FALSE)
  if (is.null(sample_size) && is.null(population_configurations) && is.null(samples) &&
      is.null(from_ts)) {
    stop("Either sample_size, samples, population_configurations or from_ts must be ",
         "specified", call. = FALSE)
  }
  samples <- .samples_factory(sample_size, samples, population_configurations)
  models <- list(.model_factory(model))
  if (!is.null(demographic_events)) {
    is_mce <- vapply(demographic_events, inherits, logical(1), "SimulationModelChange")
    model_changes <- demographic_events[is_mce]
    demographic_events <- demographic_events[!is_mce]
    current_time <- start_time %||% 0
    for (mce in model_changes) {
      last <- models[[length(models)]]
      last$duration <- mce$time - current_time
      models[[length(models)]] <- last
      models[[length(models) + 1]] <- .model_factory(mce$model)
      current_time <- mce$time
    }
  }
  is_dtwf <- inherits(models[[1]], "DiscreteTimeWrightFisher")
  is_pedigree <- any(vapply(models, inherits, logical(1), "FixedPedigree"))
  demography <- Demography$from_old_style(population_configurations,
                                          migration_matrix = migration_matrix,
                                          demographic_events = demographic_events, Ne = Ne,
                                          ignore_sample_size = TRUE)$validate()
  from_tables <- NULL
  if (is.null(from_ts)) {
    if (length(samples) < 2) stop("Sample size must be >= 2", call. = FALSE)
  } else {
    if (length(samples) > 0) stop("Cannot specify samples with from_ts", call. = FALSE)
    from_tables <- ts_tables(from_ts)
    if (demography$num_populations != length(from_tables$populations$metadata_offset) - 1) {
      stop("Mismatch in the number of populations in from_ts and simulation parameters. ",
           "The number of populations in the simulation must be equal to the number of ",
           "populations in from_ts", call. = FALSE)
    }
  }
  discrete_genome <- FALSE
  if (is.null(recombination_map)) {
    the_length <- length %||% if (is.null(from_tables)) 1 else from_tables$sequence_length
    the_rate <- recombination_rate %||% 0
    if (the_length <= 0) stop("Cannot provide non-positive sequence length", call. = FALSE)
    if (the_rate < 0) stop("Cannot provide negative recombination rate", call. = FALSE)
    recombination_map <- RateMap$uniform(the_length, the_rate)
  } else {
    if (inherits(recombination_map, "RecombinationMap")) {
      discrete_genome <- recombination_map$is_discrete
      recombination_map <- recombination_map$map
    } else if (!inherits(recombination_map, "RateMap")) {
      stop("RateMap instance required.", call. = FALSE)
    }
    if (!is.null(length) || !is.null(recombination_rate)) {
      stop("Cannot specify length/recombination_rate along with a recombination map",
           call. = FALSE)
    }
  }
  if (!is.null(from_tables) && recombination_map$sequence_length != from_tables$sequence_length) {
    stop("Recombination map and from_ts must have identical sequence_length", call. = FALSE)
  }
  if (!is.null(num_labels) && num_labels < 1) {
    stop("Must have at least one structured coalescent label", call. = FALSE)
  }
  tables <- if (is.null(from_tables)) {
    .build_initial_tables_simulate(recombination_map$sequence_length, samples, demography)
  } else {
    from_tables
  }
  additional_nodes <- NULL
  coalescing_segments_only <- TRUE
  if (record_full_arg) {
    if (is_dtwf || is_pedigree) {
      stop("Full ARG recording not supported in FixedPedigree and DTWF simulation",
           call. = FALSE)
    }
    coalescing_segments_only <- FALSE
    additional_nodes <- NodeType$RECOMBINANT + NodeType$COMMON_ANCESTOR +
      NodeType$MIGRANT + NodeType$GENE_CONVERSION
  }
  Simulator$new(
    tables = tables, recombination_map = recombination_map, models = models,
    store_migrations = record_migrations, additional_nodes = additional_nodes,
    coalescing_segments_only = coalescing_segments_only, start_time = start_time,
    end_time = end_time, num_labels = num_labels, demography = demography,
    gene_conversion_map = RateMap$uniform(recombination_map$sequence_length, 0),
    gene_conversion_tract_length = 0, discrete_genome = discrete_genome, ploidy = 2L,
    random_seed = .parse_random_seed(random_seed))
}
