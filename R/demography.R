# Demographic models (port of msprime/demography.py: Population, Demography).

.is_identifier <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x) &&
    grepl("^[\\p{L}\\p{Nl}_][\\p{L}\\p{Nl}\\p{Mn}\\p{Mc}\\p{Nd}\\p{Pc}]*$", x, perl = TRUE)
}

#' Population
#'
#' A single population in a [Demography]. Populations should be created with
#' `Demography$add_population()` rather than directly. Fields can be modified in
#' place, e.g. `demography$get("A")$initial_size <- 100`.
#'
#' @field initial_size The absolute size of the population at time zero.
#' @field growth_rate The forwards-time exponential growth rate.
#' @field name The (unique) name of the population.
#' @field description A concise description of the population.
#' @field extra_metadata A named list of extra metadata to store in the tree
#'   sequence population table.
#' @field default_sampling_time The default time at which samples are taken.
#' @field initially_active Whether the population is initially active (NULL:
#'   determined automatically).
#' @field id The zero-based ID of the population (assigned by the Demography).
#' @export
Population <- R6::R6Class("Population",
  public = list(
    initial_size = 0,
    growth_rate = 0,
    name = NULL,
    description = "",
    extra_metadata = NULL,
    default_sampling_time = NULL,
    initially_active = NULL,
    id = NULL,
    #' @description Create a population.
    #' @param initial_size,growth_rate,name,description,extra_metadata,default_sampling_time,initially_active,id
    #'   See the fields.
    initialize = function(initial_size = 0, growth_rate = 0, name = NULL,
                          description = "", extra_metadata = list(),
                          default_sampling_time = NULL, initially_active = NULL,
                          id = NULL) {
      self$initial_size <- initial_size
      self$growth_rate <- growth_rate
      self$name <- name
      self$description <- description
      self$extra_metadata <- extra_metadata %||% list()
      self$default_sampling_time <- default_sampling_time
      self$initially_active <- initially_active
      self$id <- id
    },
    #' @description Constructor arguments as a list.
    asdict = function() {
      list(initial_size = self$initial_size, growth_rate = self$growth_rate,
           name = self$name, description = self$description,
           extra_metadata = .json_object(self$extra_metadata),
           default_sampling_time = self$default_sampling_time,
           initially_active = self$initially_active, id = self$id)
    },
    #' @description Check the population is valid.
    validate = function() {
      if (!.is_number(self$initial_size)) stop("initial_size must be a number", call. = FALSE)
      if (self$initial_size < 0) stop("Negative population size", call. = FALSE)
      if (is.null(self$name)) stop("A population name must be set.", call. = FALSE)
      if (!.is_identifier(self$name)) {
        stop("A population name must be a valid identifier", call. = FALSE)
      }
      invisible(self)
    },
    #' @description Test for equality with another population.
    #' @param other Another Population.
    equals = function(other) {
      inherits(other, "Population") && isTRUE(all.equal(self$asdict(), other$asdict(),
                                                        tolerance = 0))
    },
    #' @description Print.
    #' @param ... Ignored.
    print = function(...) {
      cat(sprintf(paste0("Population(initial_size=%s, growth_rate=%s, name='%s', ",
                         "description='%s', default_sampling_time=%s, ",
                         "initially_active=%s, id=%s)\n"),
                  format(self$initial_size), format(self$growth_rate), self$name %||% "",
                  self$description %||% "", format(self$default_sampling_time %||% "NULL"),
                  format(self$initially_active %||% "NULL"), format(self$id %||% "NULL")))
      invisible(self)
    }
  )
)

.population_metadata_schema <- function() {
  list(codec = "json", type = "object",
       properties = list(name = list(type = "string"),
                         description = list(type = c("string", "null"))),
       required = c("name", "description"), additionalProperties = TRUE)
}

.encode_population_metadata <- function(population) {
  md <- list(name = population$name, description = population$description)
  extra <- population$extra_metadata
  if (length(extra) > 0) {
    clash <- intersect(names(extra), names(md))
    if (length(clash) > 0) {
      .stopf(paste0("Cannot set standard metadata key(s) [%s] using extra_metadata. ",
                    "Please set using the corresponding property of the Population ",
                    "class."), paste0("'", sort(clash), "'", collapse = ", "))
    }
    md <- c(md, extra)
  }
  charToRaw(canonical_json(md))
}

#' Demographic models
#'
#' A `Demography` defines the populations, migration rates and demographic
#' events of a model, as in Python's `msprime.Demography`. Create an empty
#' model with `Demography$new()` and add populations with `add_population()`,
#' or use one of the static constructors:
#'
#' * `Demography$isolated_model(initial_size, growth_rate = NULL)`
#' * `Demography$island_model(initial_size, migration_rate, growth_rate = NULL)`
#' * `Demography$stepping_stone_model(initial_size, migration_rate,
#'    growth_rate = NULL, boundaries = FALSE)`
#' * `Demography$from_species_tree(tree, initial_size, time_units = "gen",
#'    generation_time = NULL, growth_rate = NULL)`
#' * `Demography$from_starbeast(tree, generation_time, time_units = "myr")`
#' * `Demography$from_demes(graph)`
#' * `Demography$from_tree_sequence(ts, initial_size = 0)`
#' * `Demography$from_old_style(population_configurations, migration_matrix,
#'    demographic_events, Ne = 1, ignore_sample_size = FALSE,
#'    population_map = NULL)`
#'
#' Populations are referred to either by name or by zero-based integer ID.
#' `demography$get(id)` (or `demography[["A"]]`, `demography[[0]]`) returns
#' the [Population] object, which can be modified in place.
#'
#' @field populations The list of [Population] objects.
#' @field events The list of demographic events.
#' @field migration_matrix The backwards-time migration rate matrix; entry
#'   `[j + 1, k + 1]` is the rate at which lineages move from population `j`
#'   to population `k` (zero-based IDs).
#' @field num_populations The number of populations.
#' @field num_events The number of demographic events.
#' @export
Demography <- R6::R6Class("Demography",
  public = list(
    populations = NULL,
    events = NULL,
    migration_matrix = NULL,
    #' @description Create a demography.
    #' @param populations A list of [Population] objects.
    #' @param events A list of demographic events.
    #' @param migration_matrix The migration matrix (default all zero).
    initialize = function(populations = list(), events = list(), migration_matrix = NULL) {
      self$populations <- populations
      self$events <- events %||% list()
      N <- length(populations)
      if (is.null(migration_matrix)) {
        migration_matrix <- matrix(0, N, N)
      } else {
        migration_matrix <- as.matrix(migration_matrix)
        storage.mode(migration_matrix) <- "double"
      }
      if (!identical(dim(migration_matrix), c(N, N))) {
        stop("Migration matrix must be square and match the number of populations",
             call. = FALSE)
      }
      self$migration_matrix <- migration_matrix
      for (event in self$events) {
        if (!inherits(event, "DemographicEvent")) {
          stop("Events must be instances of DemographicEvent", call. = FALSE)
        }
      }
      for (j in seq_len(N)) {
        for (k in seq_len(j - 1)) {
          if (identical(populations[[j]], populations[[k]])) {
            stop("Population objects must be distinct", call. = FALSE)
          }
        }
      }
      for (j in seq_len(N)) {
        pop <- populations[[j]]
        if (!is.null(pop$id) && pop$id != j - 1) {
          .stopf(paste0("Population ID should be unset or set to its index within the ",
                        "Demography (expected NULL or %d, got %s)"), j - 1, pop$id)
        }
        pop$id <- j - 1
        if (is.null(pop$name)) pop$name <- sprintf("pop_%d", j - 1)
      }
      private$validate_populations()
    },
    #' @description Add a new population; its ID is the current number of
    #'   populations.
    #' @param initial_size The number of individuals at time zero.
    #' @param growth_rate The exponential growth rate (default 0).
    #' @param name The population name (default "pop_<id>").
    #' @param description A description (default "").
    #' @param extra_metadata Extra metadata (named list).
    #' @param default_sampling_time The default sampling time.
    #' @param initially_active Whether the population is initially active.
    #' @return The new [Population] (invisibly).
    add_population = function(initial_size, growth_rate = NULL, name = NULL,
                              description = NULL, extra_metadata = NULL,
                              default_sampling_time = NULL, initially_active = NULL) {
      N <- self$num_populations
      population <- Population$new(
        id = N, initial_size = initial_size, growth_rate = growth_rate %||% 0,
        name = name %||% sprintf("pop_%d", N), description = description %||% "",
        extra_metadata = extra_metadata %||% list(),
        default_sampling_time = default_sampling_time,
        initially_active = initially_active)
      # Validate before modifying, so a failed call leaves the model unchanged
      # (msprime appends the population first and leaves an invalid model).
      population$validate()
      if (population$name %in% names(self)) {
        .stopf("Duplicate population name: '%s'", population$name)
      }
      self$populations[[N + 1]] <- population
      M <- self$migration_matrix
      self$migration_matrix <- matrix(0, N + 1, N + 1)
      if (N > 0) self$migration_matrix[seq_len(N), seq_len(N)] <- M
      private$validate_populations()
      invisible(population)
    },
    #' @description Add an event object to the list of events.
    #' @param event A demographic event.
    add_event = function(event) {
      if (!inherits(event, "DemographicEvent")) {
        stop("Events must be instances of DemographicEvent", call. = FALSE)
      }
      self$events[[length(self$events) + 1]] <- event
      invisible(event)
    },
    #' @description Look up a population by name or zero-based ID.
    #' @param identifier A population name or zero-based integer ID.
    get = function(identifier) {
      if (is.character(identifier) && length(identifier) == 1) {
        for (pop in self$populations) {
          if (identical(pop$name, identifier)) return(pop)
        }
        .stopf("Population with name '%s' not found", identifier)
      }
      if (is.numeric(identifier) && length(identifier) == 1 && .isinteger(identifier)) {
        if (identifier < 0 || identifier >= self$num_populations) {
          .stopf("Population id %s out of bounds", identifier)
        }
        return(self$populations[[identifier + 1]])
      }
      stop("Keys must be either string population names or integer IDs", call. = FALSE)
    },
    #' @description Set the backwards-time migration rate from `source` to `dest`.
    #' @param source,dest Population identifiers.
    #' @param rate The per-generation migration rate.
    set_migration_rate = function(source, dest, rate) {
      source <- self$get(source)$id
      dest <- self$get(dest)$id
      if (source == dest) stop("The source and dest populations must be different",
                               call. = FALSE)
      self$migration_matrix[source + 1, dest + 1] <- rate
      invisible(self)
    },
    #' @description Set symmetric migration between all pairs in `populations`.
    #' @param populations A vector/list of population identifiers.
    #' @param rate The migration rate.
    set_symmetric_migration_rate = function(populations, rate) {
      populations <- as.list(populations)
      if (length(populations) < 2) stop("Must specify at least two populations",
                                        call. = FALSE)
      ids <- vapply(populations, function(p) self$get(p)$id, numeric(1))
      for (pair in utils::combn(ids, 2, simplify = FALSE)) {
        self$migration_matrix[pair[1] + 1, pair[2] + 1] <- rate
        self$migration_matrix[pair[2] + 1, pair[1] + 1] <- rate
      }
      invisible(self)
    },
    #' @description Add a population split: at `time`, all lineages in the
    #'   `derived` populations move to `ancestral`.
    #' @param time The time of the event.
    #' @param derived The derived population(s).
    #' @param ancestral The ancestral population.
    add_population_split = function(time, derived, ancestral) {
      derived <- as.list(derived)
      private$check_population_references(c(derived, list(ancestral)))
      pop <- self$get(ancestral)
      if (is.null(pop$initially_active)) {
        pop$initially_active <- FALSE
        if (is.null(pop$default_sampling_time)) pop$default_sampling_time <- time
      }
      self$add_event(.PopulationSplit(time, derived, ancestral))
    },
    #' @description Add an admixture event.
    #' @param time The time of the event.
    #' @param derived The admixed population.
    #' @param ancestral The ancestral populations.
    #' @param proportions The ancestry proportions (must sum to 1).
    add_admixture = function(time, derived, ancestral, proportions) {
      ancestral <- as.list(ancestral)
      private$check_population_references(c(ancestral, list(derived)))
      if (!isTRUE(abs(.psum(unlist(proportions)) - 1) <= 1e-9 * max(1, abs(.psum(unlist(proportions)))))) {
        stop("Sum of the admixture proportions must be approximately 1", call. = FALSE)
      }
      self$add_event(.Admixture(time, derived, ancestral, proportions))
    },
    #' @description Add a mass migration (pulse) event.
    #' @param time The time of the event.
    #' @param source The population lineages move from (backwards in time).
    #' @param dest The population lineages move to.
    #' @param proportion The probability that each lineage moves.
    add_mass_migration = function(time, source, dest, proportion) {
      private$check_population_references(list(source, dest))
      self$add_event(MassMigration(time, source, dest, proportion))
    },
    #' @description Change migration rates at a given time.
    #' @param time The time of the event.
    #' @param rate The new rate.
    #' @param source,dest Populations (NULL: all off-diagonal entries).
    add_migration_rate_change = function(time, rate, source = NULL, dest = NULL) {
      if (is.null(source)) {
        source <- -1
      } else {
        private$check_population_references(list(source))
      }
      if (is.null(dest)) {
        dest <- -1
      } else {
        private$check_population_references(list(dest))
      }
      self$add_event(MigrationRateChange(time = time, rate = rate, source = source,
                                         dest = dest))
    },
    #' @description Change the symmetric migration rate among populations.
    #' @param time The time of the event.
    #' @param populations Population identifiers.
    #' @param rate The new rate.
    add_symmetric_migration_rate_change = function(time, populations, rate) {
      populations <- as.list(populations)
      private$check_population_references(populations)
      self$add_event(.SymmetricMigrationRateChange(time, populations, rate))
    },
    #' @description Change population size and/or growth rate.
    #' @param time The time of the event.
    #' @param initial_size The new size (NULL: computed from previous epoch).
    #' @param growth_rate The new growth rate (NULL: unchanged).
    #' @param population The population (NULL: all populations).
    add_population_parameters_change = function(time, initial_size = NULL,
                                                growth_rate = NULL, population = NULL) {
      if (!is.null(population)) private$check_population_references(list(population))
      self$add_event(PopulationParametersChange(time, initial_size = initial_size,
                                                growth_rate = growth_rate,
                                                population = population))
    },
    #' @description Add a simple bottleneck.
    #' @param time The time of the event.
    #' @param population The population.
    #' @param proportion The probability each lineage coalesces (default 1).
    add_simple_bottleneck = function(time, population, proportion = NULL) {
      proportion <- proportion %||% 1.0
      private$check_population_references(list(population))
      self$add_event(SimpleBottleneck(time, population, proportion))
    },
    #' @description Add an instantaneous bottleneck.
    #' @param time The time of the event.
    #' @param population The population.
    #' @param strength Equivalent coalescent duration.
    add_instantaneous_bottleneck = function(time, population, strength) {
      private$check_population_references(list(population))
      self$add_event(InstantaneousBottleneck(time, population, strength))
    },
    #' @description Add a census event.
    #' @param time The time of the census.
    add_census = function(time) {
      self$add_event(CensusEvent(time))
    },
    #' @description Internal: add a population activation event.
    #' @param time The time.
    #' @param population The population.
    add_activate_population_event = function(time, population) {
      private$check_population_references(list(population))
      self$add_event(.ActivatePopulationEvent(time, population))
    },
    #' @description Sort events by time (stable).
    sort_events = function() {
      if (length(self$events) > 0) {
        times <- vapply(self$events, function(e) e$time, numeric(1))
        self$events <- self$events[order(times, method = "radix")]
      }
      invisible(self)
    },
    #' @description Validate the model and return a copy with default values
    #'   resolved.
    validate = function() {
      private$validate_populations()
      N <- self$num_populations
      if (!identical(dim(self$migration_matrix), c(as.integer(N), as.integer(N)))) {
        stop("migration matrix must be a N x N square matrix, where N is the number ",
             "of populations. The diagonal elements of this matrix must be zero.",
             call. = FALSE)
      }
      last_time <- -Inf
      for (event in self$events) {
        if (!inherits(event, "DemographicEvent")) {
          stop("Demographic events must be a list of DemographicEvent instances sorted ",
               "in non-decreasing order of time.", call. = FALSE)
        }
        if (last_time > event$time) {
          stop("Events must be time-sorted. Please use demography$sort_events() if you ",
               "add events out of order.", call. = FALSE)
        }
        last_time <- event$time
      }
      resolved <- self$deep_copy()
      for (pop in resolved$populations) {
        if (is.null(pop$default_sampling_time)) pop$default_sampling_time <- 0
        if (is.null(pop$initially_active)) pop$initially_active <- TRUE
      }
      resolved
    },
    #' @description An exact (deep) copy of this object.
    deep_copy = function() {
      pops <- lapply(self$populations, function(p) {
        np <- Population$new(initial_size = p$initial_size, growth_rate = p$growth_rate,
                             name = p$name, description = p$description,
                             default_sampling_time = p$default_sampling_time,
                             initially_active = p$initially_active, id = p$id)
        np$extra_metadata <- p$extra_metadata
        np
      })
      Demography$new(populations = pops, events = self$events,
                     migration_matrix = self$migration_matrix)
    },
    #' @description Copy of the model, optionally reordering the populations.
    #' @param populations Population identifiers giving the new order.
    copy = function(populations = NULL) {
      if (is.null(populations)) populations <- as.list(seq_len(self$num_populations) - 1)
      populations <- as.list(populations)
      if (length(populations) != self$num_populations) {
        stop("populations must have Demography.num_populations entries", call. = FALSE)
      }
      copy_demog <- Demography$new()
      for (identifier in populations) {
        pop <- self$get(identifier)
        copy_demog$add_population(
          initial_size = pop$initial_size, growth_rate = pop$growth_rate, name = pop$name,
          description = pop$description, extra_metadata = pop$extra_metadata,
          default_sampling_time = pop$default_sampling_time,
          initially_active = pop$initially_active)
      }
      N <- self$num_populations
      if (N > 1) {
        for (pair in utils::combn(seq_len(N) - 1, 2, simplify = FALSE)) {
          p1 <- copy_demog$populations[[pair[1] + 1]]
          p2 <- copy_demog$populations[[pair[2] + 1]]
          s1 <- self$get(p1$name)$id
          s2 <- self$get(p2$name)$id
          copy_demog$migration_matrix[p1$id + 1, p2$id + 1] <- self$migration_matrix[s1 + 1, s2 + 1]
          copy_demog$migration_matrix[p2$id + 1, p1$id + 1] <- self$migration_matrix[s2 + 1, s1 + 1]
        }
      }
      remap <- function(identifier) {
        if (is.numeric(identifier) && identifier == -1) return(-1)
        self$get(identifier)$name
      }
      for (event in self$events) {
        e <- event
        cls <- class(event)[1]
        if (cls %in% c("MassMigration", "MigrationRateChange")) {
          e$source <- remap(event$source)
          e$dest <- remap(event$dest)
        } else if (cls == "PopulationSplit") {
          e$derived <- lapply(event$derived, remap)
          e$ancestral <- remap(event$ancestral)
        } else if (cls == "Admixture") {
          e$ancestral <- lapply(event$ancestral, remap)
          e$derived <- remap(event$derived)
        } else if (cls %in% c("PopulationParametersChange", "SimpleBottleneck",
                              "InstantaneousBottleneck", "ActivatePopulationEvent")) {
          e$population <- remap(event$population)
        } else if (cls == "SymmetricMigrationRateChange") {
          e$populations <- lapply(event$populations, remap)
        } else if (cls != "CensusEvent") {
          stop("Event class not implemented in copy", call. = FALSE)
        }
        copy_demog$add_event(e)
      }
      copy_demog
    },
    #' @description Internal: population table (dict) for this demography.
    #' @param tables A table dict with an empty population table.
    #' @return The updated tables.
    insert_populations = function(tables) {
      stopifnot(length(tables$populations$metadata_offset) <= 1)
      md <- lapply(self$populations, .encode_population_metadata)
      packed <- .pack_raw(md)
      tables$populations <- list(metadata = packed$data, metadata_offset = packed$offset,
                                 metadata_schema = canonical_json(.population_metadata_schema()))
      tables
    },
    #' @description Internal: append populations beyond those in `tables`.
    #' @param tables A table dict.
    #' @return The updated tables.
    insert_extra_populations = function(tables) {
      schema_str <- tables$populations$metadata_schema %||% ""
      n_existing <- length(tables$populations$metadata_offset) - 1
      extra <- if (self$num_populations > n_existing) {
        self$populations[(n_existing + 1):self$num_populations]
      } else {
        list()
      }
      if (nzchar(schema_str)) {
        schema <- .json_decode(schema_str)
        if (!identical(schema$codec, "json")) {
          stop("Only JSON metadata schemas are supported for the population table",
               call. = FALSE)
        }
        properties <- schema$properties %||% list()
        additional <- schema$additionalProperties %||% TRUE
        additional <- !identical(additional, FALSE)
        name_in <- "name" %in% names(properties) || additional
        desc_in <- "description" %in% names(properties) || additional
        if (!name_in) {
          warning("The metadata schema does not have a 'name' property; population names ",
                  "will not be recorded in the output tree sequence", call. = FALSE)
        }
        if (!desc_in) {
          warning("The metadata schema does not have a 'description' property; ",
                  "population descriptions will not be recorded in the output tree ",
                  "sequence", call. = FALSE)
        }
        left_out <- character(0)
        rows <- list()
        for (population in extra) {
          md <- .json_object()
          if (name_in) md$name <- population$name
          if (desc_in) md$description <- population$description
          for (k in names(population$extra_metadata)) {
            if (k %in% names(properties) || additional) {
              md[[k]] <- population$extra_metadata[[k]]
            } else {
              left_out <- c(left_out, k)
            }
          }
          rows[[length(rows) + 1]] <- charToRaw(canonical_json(md))
        }
        tables$populations <- .append_rows(tables$populations, list(),
                                           list(metadata = rows))
        if (length(left_out) > 0) {
          warning("The metadata schema does not allow for all properties specified in ",
                  "extra_metadata: these keys have not been recorded in the output tree ",
                  "sequence: ", paste(left_out, collapse = ", "), call. = FALSE)
        }
      } else {
        warning("No metadata schema present in population table, not recording metadata",
                call. = FALSE)
        tables$populations <- .append_rows(
          tables$populations, list(), list(metadata = rep(list(raw(0)), length(extra))))
      }
      tables
    },
    #' @description List representation (used for provenance).
    asdict = function() {
      list(populations = lapply(self$populations, function(p) {
             d <- p$asdict()
             d[["__class__"]] <- "msprimeR.Population"
             d
           }),
           events = lapply(self$events, .provenance_value),
           migration_matrix = lapply(seq_len(nrow(self$migration_matrix)),
                                     function(j) I(self$migration_matrix[j, ])))
    },
    #' @description A [DemographyDebugger] for this model.
    debug = function() DemographyDebugger$new(demography = self),
    #' @description Test exact equality with another demography.
    #' @param other Another Demography.
    equals = function(other) {
      tryCatch({
        self$assert_equal(other)
        TRUE
      }, error = function(e) FALSE)
    },
    #' @description Raise an error if not exactly equal to `other`.
    #' @param other Another Demography.
    assert_equal = function(other) {
      if (!inherits(other, "Demography")) stop("Not a Demography", call. = FALSE)
      if (self$num_populations != other$num_populations) {
        stop("Number of populations differ", call. = FALSE)
      }
      for (j in seq_len(self$num_populations)) {
        if (!self$populations[[j]]$equals(other$populations[[j]])) {
          .stopf("Population %d differs", j - 1)
        }
      }
      if (!identical(unname(self$migration_matrix), unname(other$migration_matrix))) {
        stop("Migration matrices differ", call. = FALSE)
      }
      if (self$num_events != other$num_events) stop("Number of events differ", call. = FALSE)
      for (j in seq_len(self$num_events)) {
        if (!.events_equal(self$events[[j]], other$events[[j]])) {
          .stopf("Event %d differs: %s != %s", j - 1, format(self$events[[j]]),
                 format(other$events[[j]]))
        }
      }
      invisible(TRUE)
    },
    #' @description Test whether two models are equivalent up to numerical
    #'   tolerance (same epochs, population parameters, migration rates and
    #'   lineage movements).
    #' @param other Another Demography.
    #' @param rel_tol,abs_tol Tolerances.
    is_equivalent = function(other, rel_tol = NULL, abs_tol = NULL) {
      tryCatch({
        self$assert_equivalent(other, rel_tol = rel_tol, abs_tol = abs_tol)
        TRUE
      }, error = function(e) {
        if (inherits(e, "msprimeR_not_equivalent")) FALSE else stop(e)
      })
    },
    #' @description Raise an error if not equivalent to `other`.
    #' @param other Another Demography.
    #' @param rel_tol,abs_tol Tolerances.
    assert_equivalent = function(other, rel_tol = NULL, abs_tol = NULL) {
      .assert_equivalent(self, other, rel_tol %||% 1e-9, abs_tol %||% 0)
    },
    #' @description Convert to a demes graph (a resolved list structure, see
    #'   [write_demes()]).
    to_demes = function() .demography_to_demes(self),
    #' @description Print a summary of the model.
    #' @param ... Ignored.
    print = function(...) {
      resolved <- tryCatch(self$validate(), error = function(e) self)
      cat("Demography\n")
      cat(sprintf("Populations (%d):\n", resolved$num_populations))
      if (resolved$num_populations > 0) {
        df <- data.frame(
          id = vapply(resolved$populations, function(p) p$id, numeric(1)),
          name = vapply(resolved$populations, function(p) p$name, character(1)),
          description = vapply(resolved$populations, function(p) p$description %||% "",
                               character(1)),
          initial_size = vapply(resolved$populations, function(p) p$initial_size,
                                numeric(1)),
          growth_rate = vapply(resolved$populations, function(p) p$growth_rate, numeric(1)),
          default_sampling_time = vapply(resolved$populations, function(p) {
            p$default_sampling_time %||% NA_real_
          }, numeric(1)),
          stringsAsFactors = FALSE)
        print(df, row.names = FALSE)
        mm <- resolved$migration_matrix
        nms <- vapply(resolved$populations, function(p) p$name, character(1))
        dimnames(mm) <- list(nms, nms)
        if (all(mm == 0)) {
          cat("Migration matrix (all zero)\n")
        } else {
          cat("Migration matrix:\n")
          print(mm)
        }
      }
      cat(sprintf("Events (%d):\n", self$num_events))
      if (self$num_events > 0) {
        df <- data.frame(
          time = vapply(self$events, function(e) e$time, numeric(1)),
          type = vapply(self$events, function(e) .event_type_str[[class(e)[1]]],
                        character(1)),
          parameters = vapply(self$events, .event_parameters, character(1)),
          stringsAsFactors = FALSE)
        print(df, row.names = FALSE, right = FALSE)
      }
      invisible(self)
    }
  ),
  active = list(
    num_populations = function() length(self$populations),
    num_events = function() length(self$events)
  ),
  private = list(
    check_population_references = function(populations) {
      for (ref in populations) self$get(ref)
    },
    validate_populations = function() {
      names <- character(0)
      for (j in seq_along(self$populations)) {
        population <- self$populations[[j]]
        if (!inherits(population, "Population")) {
          stop("populations must be Population objects", call. = FALSE)
        }
        if (is.null(population$id) || population$id != j - 1) {
          stop("Incorrect population ID. ID values should not be updated by users. ",
               "Please use Demography$add_population to add extra populations after ",
               "initialisation.", call. = FALSE)
        }
        population$validate()
        if (population$name %in% names) {
          .stopf("Duplicate population name: '%s'", population$name)
        }
        names <- c(names, population$name)
      }
    }
  )
)

#' @export
`[[.Demography` <- function(x, i, ...) {
  if (is.character(i) && length(i) == 1) {
    pops <- .subset2(x, "populations")
    for (pop in pops) if (identical(pop$name, i)) return(pop)
    return(.subset2(x, i))
  }
  x$get(i)
}

# Supports `demography[["A"]]$initial_size <- 10`: populations are reference
# objects that are modified in place, so the replacement is a no-op for them.
#' @export
`[[<-.Demography` <- function(x, i, value) {
  if (inherits(value, "Population")) {
    if (!identical(x$get(i), value)) {
      stop("Populations cannot be replaced; modify them in place or use add_population()",
           call. = FALSE)
    }
    return(x)
  }
  if (!is.character(i) || !exists(i, envir = x, inherits = FALSE)) {
    stop("Cannot add fields to a Demography", call. = FALSE)
  }
  assign(i, value, envir = x)
  x
}

#' @export
names.Demography <- function(x) {
  vapply(.subset2(x, "populations"), function(p) p$name, character(1))
}

#' @export
length.Demography <- function(x) length(.subset2(x, "populations"))

Demography$isolated_model <- function(initial_size, growth_rate = NULL) {
  initial_size <- as.double(initial_size)
  if (!is.null(dim(initial_size)) && length(dim(initial_size)) > 1) {
    stop("The initial_size argument must a 1D array of population size values",
         call. = FALSE)
  }
  growth_rate <- if (is.null(growth_rate)) rep(0, length(initial_size)) else as.double(growth_rate)
  if (length(growth_rate) != length(initial_size)) {
    stop("If growth_rate is specified it must be a 1D array of the same length as the ",
         "population_size array", call. = FALSE)
  }
  if (any(initial_size < 0)) stop("population size values must be nonnegative.", call. = FALSE)
  if (!all(is.finite(initial_size))) stop("population size values must be finite.", call. = FALSE)
  if (!all(is.finite(growth_rate))) stop("growth_rate values must be finite.", call. = FALSE)
  pops <- lapply(seq_along(initial_size), function(j) {
    Population$new(initial_size = initial_size[j], growth_rate = growth_rate[j])
  })
  Demography$new(populations = pops)
}

Demography$island_model <- function(initial_size, migration_rate, growth_rate = NULL) {
  model <- Demography$isolated_model(initial_size, growth_rate = growth_rate)
  if (migration_rate < 0) stop("Migration rates must be non-negative", call. = FALSE)
  N <- model$num_populations
  mm <- matrix(migration_rate, N, N)
  diag(mm) <- 0
  model$migration_matrix <- mm
  model
}

Demography$stepping_stone_model <- function(initial_size, migration_rate, growth_rate = NULL,
                                            boundaries = FALSE) {
  if (!is.null(dim(initial_size)) && length(dim(initial_size)) > 1) {
    stop("Only 1D stepping stone models currently supported. Please open an issue on ",
         "GitHub if you would like 2D (or more) models", call. = FALSE)
  }
  model <- Demography$isolated_model(initial_size, growth_rate = growth_rate)
  if (migration_rate < 0) stop("Migration rates must be non-negative", call. = FALSE)
  N <- model$num_populations
  if (N > 1) {
    for (j in seq_len(N)) {
      k <- (j %% N) + 1
      model$migration_matrix[j, k] <- migration_rate
      model$migration_matrix[k, j] <- migration_rate
    }
    if (boundaries) {
      model$migration_matrix[1, N] <- 0
      model$migration_matrix[N, 1] <- 0
    }
  }
  model
}

Demography$from_tree_sequence <- function(ts, initial_size = 0) {
  tables <- ts_tables(ts)
  demography <- Demography$new()
  md <- .decode_population_metadata(tables)
  for (m in md) {
    name <- if (is.list(m)) m$name else NULL
    description <- if (is.list(m)) m$description else NULL
    demography$add_population(initial_size = initial_size, name = name,
                              description = description)
  }
  demography
}

# Decode population metadata using the table's JSON schema. Returns a list with
# one element per population: a named list if decodable, otherwise raw bytes.
.decode_population_metadata <- function(tables) {
  pt <- tables$populations
  rows <- .unpack_raw(pt$metadata, pt$metadata_offset)
  schema_str <- pt$metadata_schema %||% ""
  if (!nzchar(schema_str)) return(rows)
  schema <- tryCatch(.json_decode(schema_str), error = function(e) NULL)
  if (!identical(schema$codec, "json")) return(rows)
  lapply(rows, function(r) {
    if (length(r) == 0) return(list())
    tryCatch(.json_decode(r), error = function(e) r)
  })
}

Demography$from_species_tree <- function(tree, initial_size, time_units = "gen",
                                         generation_time = NULL, growth_rate = NULL) {
  parse_species_tree(tree, initial_size = initial_size, growth_rate = growth_rate,
                     time_units = time_units, generation_time = generation_time)
}

Demography$from_starbeast <- function(tree, generation_time, time_units = "myr") {
  parse_starbeast(tree, generation_time = generation_time, time_units = time_units)
}

Demography$from_demes <- function(graph) {
  .demography_from_demes(graph)
}

Demography$from_old_style <- function(population_configurations = NULL, migration_matrix = NULL,
                                      demographic_events = NULL, Ne = 1,
                                      ignore_sample_size = FALSE, population_map = NULL) {
  if (is.null(population_configurations)) {
    pop_configs <- list(PopulationConfiguration(initial_size = Ne))
  } else {
    pop_configs <- lapply(population_configurations, function(pc) {
      if (!inherits(pc, "PopulationConfiguration")) {
        stop("Population configurations must be a list of PopulationConfiguration ",
             "instances", call. = FALSE)
      }
      if (is.null(pc$initial_size)) pc$initial_size <- Ne
      if (!is.null(pc$sample_size) && !ignore_sample_size) {
        stop("You have specified a `sample_size` in a PopulationConfiguration object ",
             "that is to be converted into a new-style Demography object, which does ",
             "not contain any information about samples. Please use the ``samples`` ",
             "argument to sim_ancestry instead, which provides flexible options for ",
             "sampling from different populations", call. = FALSE)
      }
      pc
    })
  }
  if (is.null(population_map)) {
    .from_old_style_simple(pop_configs, migration_matrix, demographic_events)
  } else {
    .from_old_style_map_populations(pop_configs, migration_matrix, demographic_events,
                                    population_map)
  }
}

#' Population configuration (deprecated 0.x API)
#'
#' The initial configuration of a population in [simulate()]. Please use
#' [Demography] in new code.
#'
#' @param sample_size The number of sampled genomes from this population.
#' @param initial_size The population size at time zero (default: `Ne`).
#' @param growth_rate The exponential growth rate.
#' @param metadata A named list of metadata (may contain `name` and
#'   `description`).
#' @return A `PopulationConfiguration` object.
#' @export
PopulationConfiguration <- function(sample_size = NULL, initial_size = NULL,
                                    growth_rate = 0.0, metadata = NULL) {
  if (!is.null(initial_size) && initial_size < 0) {
    stop("Population size must be >= 0", call. = FALSE)
  }
  if (!is.null(sample_size) && sample_size < 0) stop("Sample size must be >= 0", call. = FALSE)
  structure(list(sample_size = sample_size, initial_size = initial_size,
                 growth_rate = growth_rate, metadata = metadata),
            class = c("PopulationConfiguration", "msprimeR_object"))
}

.add_population_from_old_style <- function(demography, pop_config, name = NULL) {
  population <- demography$add_population(name = name, initial_size = pop_config$initial_size,
                                          growth_rate = pop_config$growth_rate)
  metadata <- pop_config$metadata
  if (!is.null(metadata) && is.list(metadata)) {
    if (!is.null(metadata$name)) {
      population$name <- metadata$name
      metadata$name <- NULL
      if (!is.null(name) && name != population$name) {
        .stopf(paste0("Population name already set in old-style metadata (%s) and ",
                      "doesn't match supplied name (%s)"), name, population$name)
      }
    }
    if (!is.null(metadata$description)) {
      population$description <- metadata$description
      metadata$description <- NULL
    }
  }
  # NULL (no metadata) is distinct from an empty list for simulate().
  population$extra_metadata <- metadata
  population
}

.from_old_style_simple <- function(population_configurations, migration_matrix = NULL,
                                   demographic_events = NULL) {
  demography <- Demography$new()
  for (pc in population_configurations) {
    name <- if (is.list(pc$metadata)) pc$metadata$name else NULL
    .add_population_from_old_style(demography, pc, name)
  }
  if (!is.null(migration_matrix)) {
    mm <- if (is.list(migration_matrix)) do.call(rbind, migration_matrix) else migration_matrix
    demography$migration_matrix <- matrix(as.double(mm), nrow(as.matrix(mm)))
  }
  for (event in demographic_events) demography$add_event(event)
  demography
}

.from_old_style_map_populations <- function(population_configurations, migration_matrix,
                                            demographic_events, population_map) {
  direct_model <- .from_old_style_simple(population_configurations, migration_matrix,
                                         demographic_events)
  for (id_map in population_map) {
    vals <- unlist(id_map)
    if (anyDuplicated(vals)) stop("Population IDs in old model must be unique", call. = FALSE)
    for (old_id in vals) {
      if (old_id < 0 || old_id >= direct_model$num_populations) {
        .stopf("Bad population reference %s in old style model", old_id)
      }
    }
  }
  if (length(population_map[[1]]) != length(population_configurations)) {
    stop("The ID map for the first epoch must have entries for all populations in the ",
         "old-style model", call. = FALSE)
  }
  dbg <- suppressWarnings(direct_model$debug())
  if (dbg$num_epochs != length(population_map)) {
    .stopf(paste0("Mismatch in the number of epochs in the old style model (%d) and ",
                  "specified in the population ID map (%d)"), dbg$num_epochs,
           length(population_map))
  }
  demography <- Demography$new()
  id_map <- population_map[[1]]
  for (new_name in names(id_map)) {
    pc <- population_configurations[[id_map[[new_name]] + 1]]
    .add_population_from_old_style(demography, pc, new_name)
  }
  for (epoch_id_map in population_map[-1]) {
    for (name in names(epoch_id_map)) {
      if (!name %in% names(demography)) demography$add_population(name = name, initial_size = 0)
    }
  }
  pairs <- function(x) if (length(x) < 2) list() else utils::combn(x, 2, simplify = FALSE)
  for (pair in pairs(names(id_map))) {
    for (pk in list(pair, rev(pair))) {
      j <- id_map[[pk[1]]]
      k <- id_map[[pk[2]]]
      demography$set_migration_rate(source = pk[1], dest = pk[2],
                                    rate = direct_model$migration_matrix[j + 1, k + 1])
    }
  }
  for (epoch in dbg$epochs[-1]) {
    epoch_id_map <- population_map[[epoch$index + 1]]
    for (new_id in names(epoch_id_map)) {
      pop_state <- epoch$populations[[epoch_id_map[[new_id]] + 1]]
      demography$add_population_parameters_change(
        epoch$start_time, population = new_id, initial_size = pop_state$start_size,
        growth_rate = pop_state$growth_rate)
    }
    for (pair in pairs(names(epoch_id_map))) {
      for (pk in list(pair, rev(pair))) {
        j <- epoch_id_map[[pk[1]]]
        k <- epoch_id_map[[pk[2]]]
        demography$add_migration_rate_change(time = epoch$start_time, source = pk[1],
                                             dest = pk[2],
                                             rate = epoch$migration_matrix[j + 1, k + 1])
      }
    }
  }
  last_epoch_id_map <- population_map[[1]]
  for (epoch in dbg$epochs) {
    epoch_id_map <- population_map[[epoch$index + 1]]
    old_id_map <- stats::setNames(names(epoch_id_map), unlist(epoch_id_map))
    last_epoch_pops <- names(last_epoch_id_map)
    epoch_pops <- names(epoch_id_map)
    events <- epoch$events
    old_derived_ids <- integer(0)
    if (!setequal(last_epoch_pops, epoch_pops)) {
      ancestral <- setdiff(epoch_pops, last_epoch_pops)
      derived <- setdiff(last_epoch_pops, epoch_pops)
      if (length(ancestral) == 1) {
        demography$add_population_split(time = epoch$start_time, derived = derived,
                                        ancestral = ancestral)
        old_derived_ids <- vapply(derived, function(p) last_epoch_id_map[[p]], numeric(1))
        derived_mass_migrations <- 0
        for (event in events) {
          if (inherits(event, "MassMigration")) {
            if (event$source %in% old_derived_ids) {
              derived_mass_migrations <- derived_mass_migrations + 1
            }
            if (event$proportion != 1) {
              stop("MassMigration associated with population split with proportion != 1",
                   call. = FALSE)
            }
          }
        }
        if (derived_mass_migrations < length(old_derived_ids) - 1) {
          stop("Insufficient MassMigrations found for population split", call. = FALSE)
        }
      } else if (length(derived) == 1) {
        anc <- character(0)
        seq_props <- numeric(0)
        old_derived_id <- last_epoch_id_map[[derived]]
        for (event in events) {
          if (inherits(event, "MassMigration") && event$source == old_derived_id) {
            anc <- c(anc, old_id_map[[as.character(event$dest)]])
            seq_props <- c(seq_props, event$proportion)
          }
        }
        proportions <- .sequential_to_proportions(seq_props)
        if (isTRUE(abs(sum(proportions) - 1) <= 1e-9)) {
          if (length(anc) == 1) {
            pop <- demography$get(anc[1])
            pop$initially_active <- TRUE
            pop$default_sampling_time <- NULL
            demography$add_population_split(time = epoch$start_time, derived = derived,
                                            ancestral = anc[1])
          } else {
            demography$add_admixture(time = epoch$start_time, derived = derived,
                                     ancestral = anc, proportions = proportions)
          }
          old_derived_ids <- old_derived_id
        } else {
          stop("Admixture or single population split implied by ID map but absolute ",
               "population proportions don't sum to 1", call. = FALSE)
        }
      }
    }
    for (event in events) {
      if (inherits(event, "MassMigration")) {
        if (!event$source %in% old_derived_ids) {
          demography$add_mass_migration(time = epoch$start_time,
                                        source = old_id_map[[as.character(event$source)]],
                                        dest = old_id_map[[as.character(event$dest)]],
                                        proportion = event$proportion)
        }
      } else if (!inherits(event, c("MigrationRateChange", "PopulationParametersChange"))) {
        stop("Only MassMigration, MigrationRateChange and PopulationParametersChange ",
             "events are supported", call. = FALSE)
      }
    }
    active <- unlist(epoch_id_map)
    n <- length(epoch$populations)
    for (pair in pairs(seq_len(n) - 1)) {
      for (sd in list(pair, rev(pair))) {
        if (epoch$migration_matrix[sd[1] + 1, sd[2] + 1] != 0) {
          if (sd[1] %in% active && !sd[2] %in% active) {
            .stopf("Non zero migration from an active population (%d) to inactive (%d)",
                   sd[1], sd[2])
          }
          if (!sd[1] %in% active) {
            .warnf(paste0("Migration out of inactive population (%d) to (%d). This may ",
                          "be an error in your model."), sd[1], sd[2])
          }
        }
      }
    }
    last_epoch_id_map <- epoch_id_map
  }
  demography$sort_events()
  demography
}
