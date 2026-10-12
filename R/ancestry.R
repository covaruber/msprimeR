# Ancestry models, sample specification and the sim_ancestry() frontend
# (port of msprime/ancestry.py).

TIME_UNITS_GENERATIONS <- "generations"

#' Node flags set by msprime
#'
#' Flags recorded in the node table for nodes arising from particular events
#' (see the `additional_nodes` and `record_full_arg` arguments of
#' [sim_ancestry()]). `NodeType` is a named list of these flags which can be
#' combined with `bitwOr()`.
#'
#' @format Numeric constants.
#' @name node_flags
#' @export
NODE_IS_RE_EVENT <- 2^17
#' @rdname node_flags
#' @export
NODE_IS_CA_EVENT <- 2^18
#' @rdname node_flags
#' @export
NODE_IS_MIG_EVENT <- 2^19
#' @rdname node_flags
#' @export
NODE_IS_CEN_EVENT <- 2^20
#' @rdname node_flags
#' @export
NODE_IS_GC_EVENT <- 2^21
#' @rdname node_flags
#' @export
NODE_IS_PASS_THROUGH <- 2^22
#' @rdname node_flags
#' @export
NodeType <- list(RECOMBINANT = 2^17, COMMON_ANCESTOR = 2^18, MIGRANT = 2^19,
                 CENSUS = 2^20, GENE_CONVERSION = 2^21, PASS_THROUGH = 2^22)

#' Exit reasons of the low-level simulator
#' @format A named integer vector.
#' @export
ExitReason <- c(COALESCENCE = 0L, MAX_EVENTS = 1L, MAX_TIME = 2L, MODEL_COMPLETE = 3L)

# ---------------------------------------------------------------------------
# Ancestry models
# ---------------------------------------------------------------------------

.make_model <- function(cls, name, fields) {
  structure(c(list(duration = fields$duration), fields[setdiff(names(fields), "duration")]),
            class = c(cls, if (cls %in% c("SMCK", "BetaCoalescent", "DiracCoalescent",
                                           "SweepGenicSelection")) "ParametricAncestryModel",
                      "AncestryModel", "msprimeR_object"),
            model_name = name)
}

#' Ancestry models
#'
#' Constructors for the ancestry models available in [sim_ancestry()]. Each
#' model may be given a `duration` (in generations) when used as part of a list
#' of models.
#'
#' * `StandardCoalescent()` (`"hudson"`): the coalescent with recombination
#'   (Hudson's algorithm).
#' * `SmcApproxCoalescent()` (`"smc"`) and `SmcPrimeApproxCoalescent()`
#'   (`"smc_prime"`): legacy implementations of the SMC and SMC'.
#' * `SMCK(k)`: the general SMC(k) model.
#' * `DiscreteTimeWrightFisher()` (`"dtwf"`): discrete-time Wright-Fisher.
#' * `FixedPedigree()` (`"fixed_pedigree"`): simulation through a pedigree
#'   given in the `initial_state`.
#' * `BetaCoalescent(alpha, truncation_point)` and `DiracCoalescent(psi, c)`:
#'   multiple merger coalescents.
#' * `SweepGenicSelection(position, start_frequency, end_frequency, s, dt)`:
#'   a selective sweep.
#'
#' @param duration The maximum time (in generations) for which this model runs.
#' @param k The SMC(k) hull offset.
#' @param alpha,truncation_point Beta coalescent parameters.
#' @param psi,c Dirac coalescent parameters.
#' @param position,start_frequency,end_frequency,s,dt Sweep parameters.
#' @return An `AncestryModel` object.
#' @name ancestry_models
NULL

#' @rdname ancestry_models
#' @export
StandardCoalescent <- function(duration = NULL) {
  .make_model("StandardCoalescent", "hudson", list(duration = duration))
}

#' @rdname ancestry_models
#' @export
SmcApproxCoalescent <- function(duration = NULL) {
  .make_model("SmcApproxCoalescent", "smc", list(duration = duration))
}

#' @rdname ancestry_models
#' @export
SmcPrimeApproxCoalescent <- function(duration = NULL) {
  .make_model("SmcPrimeApproxCoalescent", "smc_prime", list(duration = duration))
}

#' @rdname ancestry_models
#' @export
SMCK <- function(k, duration = NULL) {
  .make_model("SMCK", "smc_k", list(duration = duration, k = k))
}

#' @rdname ancestry_models
#' @export
DiscreteTimeWrightFisher <- function(duration = NULL) {
  .make_model("DiscreteTimeWrightFisher", "dtwf", list(duration = duration))
}

#' @rdname ancestry_models
#' @export
FixedPedigree <- function(duration = NULL) {
  .make_model("FixedPedigree", "fixed_pedigree", list(duration = duration))
}

#' @rdname ancestry_models
#' @export
BetaCoalescent <- function(alpha = NULL, truncation_point = .Machine$double.xmax,
                           duration = NULL) {
  .make_model("BetaCoalescent", "beta", list(duration = duration, alpha = alpha,
                                              truncation_point = truncation_point))
}

#' @rdname ancestry_models
#' @export
DiracCoalescent <- function(psi = NULL, c = NULL, duration = NULL) {
  .make_model("DiracCoalescent", "dirac", list(duration = duration, psi = psi, c = c))
}

#' @rdname ancestry_models
#' @export
SweepGenicSelection <- function(position = NULL, start_frequency = NULL,
                                end_frequency = NULL, s = NULL, dt = NULL,
                                duration = NULL) {
  .make_model("SweepGenicSelection", "sweep_genic_selection",
              list(duration = duration, position = position,
                   start_frequency = start_frequency, end_frequency = end_frequency,
                   s = s, dt = dt))
}

#' @export
print.AncestryModel <- function(x, ...) {
  params <- unclass(x)
  attr(params, "model_name") <- NULL
  params <- params[!vapply(params, is.null, logical(1))]
  s <- paste(names(params), vapply(params, format, character(1)), sep = "=", collapse = ", ")
  cat(sprintf("%s(%s)\n", class(x)[1], s))
  invisible(x)
}

.model_as_lowlevel <- function(model) {
  d <- unclass(model)
  attr(d, "model_name") <- NULL
  d$duration <- NULL
  d <- d[!vapply(d, is.null, logical(1))]
  c(list(name = attr(model, "model_name")), d)
}

.model_factory <- function(model) {
  model_map <- list(hudson = StandardCoalescent, smc = SmcApproxCoalescent,
                    smc_prime = SmcPrimeApproxCoalescent, dtwf = DiscreteTimeWrightFisher,
                    fixed_pedigree = FixedPedigree)
  if (is.null(model)) return(StandardCoalescent())
  if (is.character(model) && length(model) == 1) {
    lower <- tolower(model)
    if (!lower %in% names(model_map)) {
      .stopf("Model '%s' unknown. Choose from %s", model,
             paste0("['", paste(names(model_map), collapse = "', '"), "']"))
    }
    return(model_map[[lower]]())
  }
  if (!inherits(model, "AncestryModel")) {
    stop("Ancestry model must be a string or an instance of AncestryModel", call. = FALSE)
  }
  model
}

.parse_model_arg <- function(model_arg) {
  if (is.list(model_arg) && !inherits(model_arg, "AncestryModel")) {
    if (length(model_arg) < 1) stop("Must specify at least one AncestryModel", call. = FALSE)
    return(lapply(model_arg, .model_factory))
  }
  if (is.character(model_arg) && length(model_arg) > 1) {
    return(lapply(model_arg, .model_factory))
  }
  list(.model_factory(model_arg))
}

#' Simulation model change (deprecated)
#'
#' Demographic-event style specification of a change of ancestry model, for
#' use with [simulate()]. In new code pass a list of models to
#' [sim_ancestry()].
#'
#' @param time The time of the model change.
#' @param model The new model (string or `AncestryModel`).
#' @return A `SimulationModelChange` object.
#' @export
SimulationModelChange <- function(time = NULL, model = NULL) {
  structure(list(time = time, model = model),
            class = c("SimulationModelChange", "msprimeR_object"))
}

#' Sample sets
#'
#' Specifies a set of exchangeable sample individuals for [sim_ancestry()].
#'
#' @param num_samples The number of k-ploid sample individuals.
#' @param population The population (name or zero-based ID).
#' @param time The sampling time (default: the population's default sampling
#'   time).
#' @param ploidy The ploidy of the sampled individuals (default: the
#'   simulation ploidy).
#' @return A `SampleSet` object.
#' @export
SampleSet <- function(num_samples, population = NULL, time = NULL, ploidy = NULL) {
  structure(list(num_samples = num_samples, population = population, time = time,
                 ploidy = ploidy), class = c("SampleSet", "msprimeR_object"))
}

#' @export
print.SampleSet <- function(x, ...) {
  f <- function(v) if (is.null(v)) "NULL" else format(v)
  cat(sprintf("SampleSet(num_samples=%s, population=%s, time=%s, ploidy=%s)\n",
              f(x$num_samples), f(x$population), f(x$time), f(x$ploidy)))
  invisible(x)
}

.Sample <- function(population, time) {
  structure(list(population = population, time = time), class = c("Sample", "msprimeR_object"))
}

#' Legacy sample specification
#'
#' A sample genome in a given population at a given time, for [simulate()].
#'
#' @param population The zero-based population ID.
#' @param time The sampling time.
#' @return A `Sample` object.
#' @export
Sample <- function(population, time) .Sample(population, time)

# ---------------------------------------------------------------------------
# Low-level simulator wrapper
# ---------------------------------------------------------------------------

Simulator <- R6::R6Class("Simulator",
  cloneable = FALSE,
  public = list(
    xp = NULL,
    end_time = Inf,
    models = NULL,
    demography = NULL,
    gene_conversion_map = NULL,
    missing_intervals = NULL,
    initialize = function(tables, recombination_map, gene_conversion_map,
                          gene_conversion_tract_length, discrete_genome, ploidy,
                          demography, random_seed, models = NULL,
                          store_migrations = FALSE, additional_nodes = NULL,
                          coalescing_segments_only = TRUE, start_time = NULL,
                          end_time = NULL, num_labels = NULL, stop_at_local_mrca = TRUE) {
      num_samples <- length(tables$nodes$flags)
      block_size <- 64 * 1024
      if (is.null(num_labels)) num_labels <- private$choose_num_labels(models)
      additional_nodes <- additional_nodes %||% 0
      pop_config <- lapply(demography$populations, function(p) {
        list(initial_size = p$initial_size, growth_rate = p$growth_rate,
             initially_active = p$initially_active)
      })
      events <- lapply(demography$events, .ll_event, demography = demography)
      ll_recomb_map <- private$resolve_missing_intervals(recombination_map)
      if (gene_conversion_map$num_intervals != 1) {
        stop("Only uniform gene conversion maps are currently supported", call. = FALSE)
      }
      params <- list(
        avl_node_block_size = block_size,
        segment_block_size = max(block_size, num_samples),
        node_mapping_block_size = block_size,
        store_migrations = store_migrations,
        start_time = if (is.null(start_time)) -1 else start_time,
        additional_nodes = additional_nodes,
        coalescing_segments_only = coalescing_segments_only,
        num_labels = num_labels,
        gene_conversion_rate = gene_conversion_map$rate[1],
        gene_conversion_tract_length = gene_conversion_tract_length,
        discrete_genome = discrete_genome,
        ploidy = ploidy,
        stop_at_local_mrca = stop_at_local_mrca)
      self$xp <- .Call(msprimer_sim_new, tables, random_seed, ll_recomb_map, pop_config,
                       demography$migration_matrix, events, NULL, params)
      self$end_time <- if (is.null(end_time)) Inf else end_time
      self$models <- models
      self$demography <- demography
      self$gene_conversion_map <- gene_conversion_map
    },
    state = function() .Call(msprimer_sim_state, self$xp),
    free = function() invisible(.Call(msprimer_sim_free, self$xp)),
    time = function() .Call(msprimer_sim_state, self$xp)$time,
    set_model = function(model) .Call(msprimer_sim_set_model, self$xp, model),
    get_model = function() .Call(msprimer_sim_get_model, self$xp),
    debug_demography = function() .Call(msprimer_sim_debug_demography, self$xp),
    compute_population_size = function(population, time) {
      .Call(msprimer_sim_compute_population_size, self$xp, population, time)
    },
    tables = function() .Call(msprimer_sim_get_tables, self$xp),
    finalise_tables = function() .Call(msprimer_sim_finalise_tables, self$xp),
    reset = function() .Call(msprimer_sim_reset, self$xp),
    ancestors = function() .Call(msprimer_sim_ancestors, self$xp),
    breakpoints = function() .Call(msprimer_sim_breakpoints, self$xp),
    print_state = function() {
      f <- tempfile()
      on.exit(unlink(f))
      .Call(msprimer_sim_print_state, self$xp, f)
      paste(readLines(f), collapse = "\n")
    },
    run_until = function(end_time, event_chunk = NULL, debug_func = NULL) {
      event_chunk <- event_chunk %||% 10^4
      if (event_chunk <= 0) stop("Must have at least 1 event per chunk", call. = FALSE)
      ret <- ExitReason[["MAX_EVENTS"]]
      while (ret == ExitReason[["MAX_EVENTS"]]) {
        ret <- .Call(msprimer_sim_run, self$xp, end_time, event_chunk)
        if (self$time() > end_time) {
          .stopf("Model %s does not support interruption. Please open an issue on GitHub",
                 self$get_model()$name)
        }
        if (!is.null(debug_func)) debug_func(self)
      }
      ret
    },
    run = function(event_chunk = NULL, debug_func = NULL) {
      for (model in self$models) {
        self$set_model(.model_as_lowlevel(model))
        duration <- if (is.null(model$duration)) Inf else model$duration
        if (duration < 0) stop("Model durations must be >= 0", call. = FALSE)
        end_time <- min(self$time() + duration, self$end_time)
        exit_reason <- self$run_until(end_time, event_chunk, debug_func)
        if (exit_reason == ExitReason[["COALESCENCE"]] || self$time() == self$end_time) break
      }
      self$finalise_tables()
      invisible(self)
    },
    run_replicates = function(num_replicates, mutation_rate = NULL,
                              provenance_dict = NULL, only_last = FALSE) {
      encoded <- NULL
      placeholder <- "@@_REPLICATE_INDEX_@@"
      if (!is.null(provenance_dict)) {
        provenance_dict$parameters$replicate_index <- placeholder
        encoded <- .json_encode_provenance(provenance_dict)
      }
      results <- vector("list", if (only_last) min(1, num_replicates) else num_replicates)
      for (j in seq_len(num_replicates)) {
        self$run()
        if (!is.null(mutation_rate)) {
          # As in msprime, these mutations use the simulator's own RNG and tables.
          .Call(msprimer_sim_simple_mutate, self$xp, mutation_rate)
        }
        if (!only_last || j == num_replicates) {
          tables <- self$tables()
          if (nrow(self$missing_intervals) > 0) {
            tables <- .Call(msprimer_tables_delete_intervals, tables,
                            as.double(self$missing_intervals[, 1]),
                            as.double(self$missing_intervals[, 2]))
          }
          if (!is.null(encoded)) {
            rec <- sub(paste0("\"", placeholder, "\""), as.character(j - 1), encoded,
                       fixed = TRUE)
            tables <- .add_provenance(tables, rec)
          }
          results[[if (only_last) 1 else j]] <- .make_result(tables)
        }
        self$reset()
      }
      results
    }
  ),
  private = list(
    choose_num_labels = function(models) {
      num_labels <- 1
      for (model in models) if (inherits(model, "SweepGenicSelection")) num_labels <- 2
      num_labels
    },
    resolve_missing_intervals = function(recombination_map) {
      self$missing_intervals <- recombination_map$missing_intervals()
      msg <- paste0("Missing regions of the genome other than the flanks are currently ",
                    "not supported. Please see https://github.com/tskit-dev/msprime/issues/1604")
      mi <- self$missing_intervals
      if (nrow(mi) > 2) stop(msg, call. = FALSE)
      for (j in seq_len(nrow(mi))) {
        if (!(mi[j, 1] == 0 || mi[j, 2] == recombination_map$sequence_length)) {
          stop(msg, call. = FALSE)
        }
      }
      rate <- recombination_map$rate
      rate[recombination_map$missing] <- 0
      list(position = recombination_map$position, rate = rate)
    }
  )
)

# ---------------------------------------------------------------------------
# sim_ancestry
# ---------------------------------------------------------------------------

.insert_sample_sets <- function(sample_sets, demography, default_ploidy, tables) {
  for (ss in sample_sets) {
    n <- ss$num_samples
    population <- demography$get(ss$population)
    time <- ss$time %||% population$default_sampling_time
    ploidy <- ss$ploidy %||% default_ploidy
    n_ind <- length(tables$individuals$flags)
    node_individual <- n_ind + rep(seq_len(n) - 1L, each = ploidy)
    tables$individuals <- .append_rows(tables$individuals, list(flags = rep(0, n)))
    N <- n * ploidy
    tables$nodes <- .append_rows(tables$nodes, list(
      flags = rep(1, N), time = rep(as.double(time), N),
      population = rep(as.integer(population$id), N),
      individual = as.integer(node_individual)))
  }
  tables
}

.parse_sample_sets <- function(sample_sets, demography) {
  sample_sets <- lapply(sample_sets, function(ss) {
    if (!inherits(ss, "SampleSet")) stop("msprimeR SampleSet object required", call. = FALSE)
    if (!.isinteger(ss$num_samples)) {
      stop("The number of samples to draw from a population must be an integer",
           call. = FALSE)
    }
    ss$num_samples <- as.integer(ss$num_samples)
    if (ss$num_samples < 0) stop("Number of samples cannot be negative", call. = FALSE)
    if (is.null(ss$population)) {
      if (demography$num_populations == 1) {
        ss$population <- 0
      } else {
        stop("Must specify a SampleSet population in multipopulation models", call. = FALSE)
      }
    }
    ss
  })
  if (sum(vapply(sample_sets, function(ss) ss$num_samples, numeric(1))) == 0) {
    stop("Zero samples specified", call. = FALSE)
  }
  sample_sets
}

.parse_samples <- function(samples, demography) {
  if (inherits(samples, "SampleSet")) {
    sample_sets <- list(samples)
  } else if (is.list(samples) && is.null(names(samples))) {
    sample_sets <- samples
  } else if ((is.list(samples) || is.numeric(samples)) && !is.null(names(samples))) {
    sample_sets <- lapply(names(samples), function(key) {
      pop <- if (key %in% names(demography)) key else {
        id <- suppressWarnings(as.numeric(key))
        if (is.na(id)) key else id
      }
      SampleSet(samples[[key]], population = pop)
    })
  } else if (.isinteger(samples)) {
    if (demography$num_populations != 1) {
      stop("Numeric samples can only be used in single population models. Please use ",
           "either the mapping form c(pop = num_samples) or a list of SampleSet objects.",
           call. = FALSE)
    }
    sample_sets <- list(SampleSet(samples))
  } else {
    stop("The value of 'samples' cannot be interpreted as sample specification. ",
         "Samples must either be a single integer, a named vector/list that maps ",
         "populations to the number of samples for that population, or a list of ",
         "SampleSet objects.", call. = FALSE)
  }
  .parse_sample_sets(sample_sets, demography)
}

.flag_set <- function(value, flag) bitwAnd(as.integer(value / 2^17), as.integer(flag / 2^17)) > 0

.parse_sim_ancestry <- function(samples = NULL, sequence_length = NULL,
                                recombination_rate = NULL, gene_conversion_rate = NULL,
                                gene_conversion_tract_length = NULL, discrete_genome = NULL,
                                population_size = NULL, demography = NULL, ploidy = NULL,
                                model = NULL, initial_state = NULL, start_time = NULL,
                                end_time = NULL, record_migrations = NULL,
                                record_full_arg = NULL, additional_nodes = NULL,
                                coalescing_segments_only = NULL, num_labels = NULL,
                                random_seed = NULL, stop_at_local_mrca = NULL,
                                init_for_debugger = FALSE) {
  start_time <- if (is.null(start_time)) 0 else .check_number(start_time, "start_time")
  end_time <- if (is.null(end_time)) Inf else .check_number(end_time, "end_time")
  discrete_genome <- .parse_flag(discrete_genome, default = TRUE)
  record_full_arg <- .parse_flag(record_full_arg, default = FALSE)

  models <- .parse_model_arg(model)
  is_dtwf <- inherits(models[[1]], "DiscreteTimeWrightFisher")
  is_pedigree <- any(vapply(models, inherits, logical(1), "FixedPedigree"))
  is_smck <- any(vapply(models, inherits, logical(1), "SMCK"))

  if (is_smck && !is.null(gene_conversion_rate)) {
    stop("Gene conversion is not supported for the SMCK model. Please refer to issue ",
         "#2399 on GitHub for details.", call. = FALSE)
  }
  if (record_full_arg) {
    if (!is.null(coalescing_segments_only)) {
      stop("Cannot set both coalescing_segments_only and record_full_arg.", call. = FALSE)
    }
    if (!is.null(additional_nodes)) {
      stop("Cannot set both record_full_arg and additional_nodes.", call. = FALSE)
    }
    coalescing_segments_only <- FALSE
    if (is_dtwf || is_pedigree) {
      stop("Full ARG recording not supported in FixedPedigree and DTWF simulation",
           call. = FALSE)
    }
    additional_nodes <- NodeType$RECOMBINANT + NodeType$COMMON_ANCESTOR +
      NodeType$MIGRANT + NodeType$GENE_CONVERSION
  }
  coalescing_segments_only <- .parse_flag(coalescing_segments_only, default = TRUE)
  if (is.null(additional_nodes)) {
    additional_nodes <- 0
  } else {
    if (!.is_number(additional_nodes) || additional_nodes < 0 ||
        additional_nodes %% 2^17 != 0 || additional_nodes >= 2^23) {
      stop("additional_nodes should be a combination of msprimeR NodeType flags.",
           call. = FALSE)
    }
    if (additional_nodes > 0 && coalescing_segments_only) {
      stop("Set coalescing_segments_only to FALSE when recording additional_nodes.",
           call. = FALSE)
    }
    if (.flag_set(additional_nodes, NodeType$MIGRANT)) {
      if (is_dtwf) {
        stop("Recording MIGRANT nodes is currently not supported in DTWF simulation.",
             call. = FALSE)
      }
      if (is_pedigree) {
        stop("Recording MIGRANT nodes is not supported in FixedPedigree simulation.",
             call. = FALSE)
      }
      if (any(vapply(models, inherits, logical(1), "SweepGenicSelection"))) {
        stop("Recording MIGRANT nodes is not supported in SweepGenicSelection simulation.",
             call. = FALSE)
      }
    }
  }
  record_migrations <- .parse_flag(record_migrations, default = FALSE)

  if (!is.null(initial_state)) {
    initial_state <- ts_tables(initial_state)
    if (identical(initial_state$time_units, "uncalibrated")) {
      stop("Cannot use a tree sequence with uncalibrated time_units as the initial state",
           call. = FALSE)
    }
    if (!identical(initial_state$time_units, TIME_UNITS_GENERATIONS)) {
      warning(structure(class = c("TimeUnitsMismatchWarning", "warning", "condition"),
        list(message = paste0(
          "The initial_state has time_units=", initial_state$time_units,
          " but time is measured in generations in msprime. This may lead to ",
          "significant discrepancies between the timescales. If you wish to suppress ",
          "this warning, use suppressWarnings(..., classes = 'TimeUnitsMismatchWarning')"),
          call = NULL)))
    }
  }

  if (is.null(sequence_length)) {
    if (!is.null(initial_state)) {
      sequence_length <- initial_state$sequence_length
    } else if (is.null(recombination_rate) && is.null(gene_conversion_rate)) {
      sequence_length <- 1
    } else if (inherits(recombination_rate, "RateMap")) {
      sequence_length <- recombination_rate$sequence_length
    } else if (inherits(gene_conversion_rate, "RateMap")) {
      sequence_length <- gene_conversion_rate$sequence_length
    } else {
      stop("A sequence_length value must be specified. This can be either via the ",
           "the sequence_length parameter itself, of implicitly through using a RateMap ",
           "instance for the recombination_rate or gene_conversion_rate parameters, or ",
           "via the initial_state tables. ", call. = FALSE)
    }
  } else {
    sequence_length <- .check_number(sequence_length, "sequence_length")
  }
  if (discrete_genome && floor(sequence_length) != sequence_length) {
    stop("Must have integer sequence length with discrete_genome=TRUE", call. = FALSE)
  }
  recombination_map <- .parse_rate_map(recombination_rate, sequence_length, "recombination")
  gene_conversion_map <- .parse_rate_map(gene_conversion_rate, sequence_length,
                                         "gene conversion")
  if (is.null(gene_conversion_tract_length)) {
    if (is.null(gene_conversion_rate)) {
      gene_conversion_tract_length <- 1
    } else {
      stop("Must specify tract length when simulating gene conversion", call. = FALSE)
    }
  } else {
    if (is.null(gene_conversion_rate)) {
      stop("Must specify gene conversion rate along with tract length", call. = FALSE)
    }
    gene_conversion_tract_length <- .check_number(gene_conversion_tract_length,
                                                  "gene_conversion_tract_length")
  }
  ploidy <- ploidy %||% 2
  if (!.isinteger(ploidy)) stop("ploidy must be an integer", call. = FALSE)
  ploidy <- as.integer(ploidy)
  if (ploidy < 1) stop("ploidy must be >= 1", call. = FALSE)

  if (is_pedigree) {
    if (!is.null(demography)) {
      stop("Cannot specify demography for FixedPedigree simulation", call. = FALSE)
    }
    if (!is.null(population_size)) {
      stop("Cannot specify population_size for FixedPedigree simulation", call. = FALSE)
    }
    if (is.null(initial_state)) {
      stop("Must specify an input pedigree using the ``initial_state`` argument for a ",
           "FixedPedigree simulation", call. = FALSE)
    }
    if (ploidy != 2) stop("Fixed pedigree simulations must have ploidy=2", call. = FALSE)
    if (gene_conversion_map$total_mass != 0) {
      stop("Gene conversion not supported in FixedPedigree simulation", call. = FALSE)
    }
    if (record_migrations) {
      stop("Migration recording not supported in FixedPedigree simulation", call. = FALSE)
    }
    if (start_time > 0) {
      stop("Cannot specify start_time in FixedPedigree simulation", call. = FALSE)
    }
    if (length(models) > 1) {
      stop("Cannot use FixedPedigree simulation in conjunction with other models",
           call. = FALSE)
    }
    demography <- Demography$from_tree_sequence(initial_state, initial_size = 1)
  }

  if (is.null(demography)) {
    if (is_dtwf && is.null(population_size)) {
      stop("When using the DTWF model, the population size must be set explicitly, ",
           "either using the population_size or demography arguments.", call. = FALSE)
    }
    if (!is.null(initial_state)) {
      if (is.null(population_size)) {
        stop("Must specify either a demography object or a population_size (for single ",
             "population models) when providing an initial_state.", call. = FALSE)
      }
      if (length(initial_state$populations$metadata_offset) - 1 > 1) {
        stop("Must specify demography for initial_state with > 1 population", call. = FALSE)
      }
    }
    population_size <- if (is.null(population_size)) 1 else
      .check_number(population_size, "population_size")
    demography <- Demography$isolated_model(population_size)
  } else if (inherits(demography, "Demography")) {
    if (!is.null(population_size)) {
      stop("Cannot specify demography and population size", call. = FALSE)
    }
  } else {
    stop("demography argument must be an instance of msprimeR Demography", call. = FALSE)
  }
  demography <- demography$validate()

  if (is.null(initial_state)) {
    if (is.null(samples) && !init_for_debugger) {
      stop("Either the samples or initial_state arguments must be provided", call. = FALSE)
    }
    initial_state <- .empty_tables(sequence_length)
    initial_state$time_units <- TIME_UNITS_GENERATIONS
    initial_state <- demography$insert_populations(initial_state)
    if (!init_for_debugger) {
      sample_sets <- .parse_samples(samples, demography)
      initial_state <- .insert_sample_sets(sample_sets, demography, ploidy, initial_state)
    }
  } else {
    if (!is.null(samples)) stop("Cannot specify both samples and initial_state", call. = FALSE)
    if (sequence_length != initial_state$sequence_length) {
      stop("The initial_state sequence length must be consistent with the value derived ",
           "from either the sequence_length, recombination_rate or gene_conversion_rate ",
           "parameters.", call. = FALSE)
    }
    num_pops <- length(initial_state$populations$metadata_offset) - 1
    if (num_pops == 0) {
      stop("initial_state tables must define at least one population.", call. = FALSE)
    }
    demography_check <- Demography$from_tree_sequence(initial_state)
    if (demography$num_populations < demography_check$num_populations) {
      .stopf(paste0("Input demography must have at least as many populations as the ",
                    "initial state population table: %d < %d"),
             demography$num_populations, demography_check$num_populations)
    }
    for (j in seq_len(demography_check$num_populations)) {
      n1 <- demography$populations[[j]]$name
      n2 <- demography_check$populations[[j]]$name
      if (n1 != n2) {
        .stopf(paste0("Population names in the input demography and the initial state ",
                      "population table must be equal: %s \u2260 %s"), n1, n2)
      }
    }
    initial_state <- demography$insert_extra_populations(initial_state)
  }
  random_seed <- .parse_random_seed(random_seed)
  stop_at_local_mrca <- stop_at_local_mrca %||% TRUE

  Simulator$new(
    tables = initial_state, recombination_map = recombination_map,
    gene_conversion_map = gene_conversion_map,
    gene_conversion_tract_length = gene_conversion_tract_length,
    discrete_genome = discrete_genome, ploidy = ploidy, demography = demography,
    models = models, store_migrations = record_migrations,
    additional_nodes = additional_nodes, coalescing_segments_only = coalescing_segments_only,
    start_time = start_time, end_time = end_time, num_labels = num_labels,
    random_seed = random_seed, stop_at_local_mrca = stop_at_local_mrca)
}

.parse_replicate_index <- function(replicate_index, random_seed, num_replicates) {
  if (is.null(replicate_index)) return(NULL)
  if (is.null(random_seed)) {
    stop("Cannot specify the replicate_index without a random_seed", call. = FALSE)
  }
  if (!is.null(num_replicates)) {
    stop("Cannot specify the replicate_index as well as num_replicates", call. = FALSE)
  }
  replicate_index <- as.integer(replicate_index)
  if (replicate_index < 0) stop("Cannot specify negative replicate_index.", call. = FALSE)
  replicate_index
}

.wrap_replicates <- function(simulator, num_replicates, replicate_index, provenance_dict,
                             mutation_rate = NULL) {
  if (is.null(num_replicates) && is.null(replicate_index)) {
    replicate_index <- 0L
  }
  if (!is.null(replicate_index)) {
    # Earlier replicates must still be simulated to reach the requested one.
    res <- simulator$run_replicates(replicate_index + 1, mutation_rate = mutation_rate,
                                    provenance_dict = provenance_dict, only_last = TRUE)
    return(res[[1]])
  }
  if (!.isinteger(num_replicates) || num_replicates < 0) {
    stop("num_replicates must be a non-negative integer", call. = FALSE)
  }
  simulator$run_replicates(num_replicates, mutation_rate = mutation_rate,
                           provenance_dict = provenance_dict)
}

.parameters_dict <- function(command, args) {
  c(list(command = command), args)
}

#' Simulate ancestral histories
#'
#' Simulates an ancestral process described by the specified model, demography
#' and samples, returning a tree sequence (or a list of replicate tree
#' sequences). This is the R equivalent of Python's `msprime.sim_ancestry()`
#' and is driven by the same C library: given the same arguments and
#' `random_seed`, it produces identical results.
#'
#' @param samples The sampled individuals: a single integer (number of
#'   individuals in a single population model), a named vector or list mapping
#'   populations (names or zero-based IDs) to numbers of individuals, or a list
#'   of [SampleSet()] objects.
#' @param demography A [Demography] (default: a single population of size
#'   `population_size`).
#' @param sequence_length The length of the genome.
#' @param discrete_genome If TRUE (default), breakpoints occur at integer
#'   positions.
#' @param recombination_rate A single rate or a [RateMap].
#' @param gene_conversion_rate The gene conversion rate (requires
#'   `gene_conversion_tract_length`).
#' @param gene_conversion_tract_length The mean gene conversion tract length.
#' @param population_size The size of the default single population.
#' @param ploidy The number of genomes per sampled individual (default 2).
#' @param model An ancestry model, model name, or list of models (see
#'   [ancestry_models]).
#' @param initial_state A tree sequence (or tables) to start the simulation
#'   from.
#' @param start_time,end_time Start and end times of the simulation.
#' @param record_migrations If TRUE, record migrations in the migration table.
#' @param record_full_arg If TRUE, record all nodes of the ARG.
#' @param additional_nodes Node types (see [NodeType]) for which all ancestry is
#'   retained, combined with `+` or `bitwOr()`.
#' @param coalescing_segments_only If FALSE, retain ancestry for coalescent
#'   nodes across the full sequence.
#' @param num_labels Number of structured coalescent labels (internal).
#' @param random_seed The random seed (an integer in `[1, 2^32 - 1]`). If NULL,
#'   a seed is drawn from R's random number generator, so `set.seed()` makes
#'   simulations reproducible.
#' @param num_replicates If specified, return a list of this many replicates.
#' @param replicate_index Return only the replicate with this (zero-based)
#'   index, which must be used along with `random_seed`.
#' @param record_provenance If TRUE (default), record a provenance entry.
#' @param stop_at_local_mrca If TRUE (default), stop simulating each tree at its
#'   local MRCA.
#' @return A tree sequence (see [msprime_ts]), or a list of them if
#'   `num_replicates` is given.
#' @examples
#' ts <- sim_ancestry(samples = 5, population_size = 1e4, sequence_length = 1e5,
#'                    recombination_rate = 1e-8, random_seed = 42)
#' @export
sim_ancestry <- function(samples = NULL, demography = NULL, sequence_length = NULL,
                         discrete_genome = NULL, recombination_rate = NULL,
                         gene_conversion_rate = NULL, gene_conversion_tract_length = NULL,
                         population_size = NULL, ploidy = NULL, model = NULL,
                         initial_state = NULL, start_time = NULL, end_time = NULL,
                         record_migrations = NULL, record_full_arg = NULL,
                         additional_nodes = NULL, coalescing_segments_only = NULL,
                         num_labels = NULL, random_seed = NULL, num_replicates = NULL,
                         replicate_index = NULL, record_provenance = NULL,
                         stop_at_local_mrca = NULL) {
  record_provenance <- record_provenance %||% TRUE
  replicate_index <- .parse_replicate_index(replicate_index, random_seed, num_replicates)
  random_seed <- .parse_random_seed(random_seed)
  provenance_dict <- NULL
  if (record_provenance) {
    provenance_dict <- .get_provenance_dict(.parameters_dict("sim_ancestry", list(
      samples = samples, demography = demography, sequence_length = sequence_length,
      discrete_genome = discrete_genome, recombination_rate = recombination_rate,
      gene_conversion_rate = gene_conversion_rate,
      gene_conversion_tract_length = gene_conversion_tract_length,
      population_size = population_size, ploidy = ploidy, model = model,
      initial_state = initial_state, start_time = start_time, end_time = end_time,
      record_migrations = record_migrations, record_full_arg = record_full_arg,
      additional_nodes = additional_nodes,
      coalescing_segments_only = coalescing_segments_only, num_labels = num_labels,
      random_seed = random_seed, stop_at_local_mrca = stop_at_local_mrca)))
  }
  sim <- .parse_sim_ancestry(
    samples = samples, sequence_length = sequence_length,
    recombination_rate = recombination_rate, gene_conversion_rate = gene_conversion_rate,
    gene_conversion_tract_length = gene_conversion_tract_length,
    discrete_genome = discrete_genome, population_size = population_size,
    demography = demography, ploidy = ploidy, model = model, initial_state = initial_state,
    start_time = start_time, end_time = end_time, record_migrations = record_migrations,
    record_full_arg = record_full_arg, additional_nodes = additional_nodes,
    coalescing_segments_only = coalescing_segments_only, num_labels = num_labels,
    random_seed = random_seed, stop_at_local_mrca = stop_at_local_mrca)
  on.exit(sim$free(), add = TRUE)
  .wrap_replicates(sim, num_replicates, replicate_index, provenance_dict)
}
