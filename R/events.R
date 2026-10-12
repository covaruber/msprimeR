# Demographic events (port of the event classes in msprime/demography.py).
#
# Events are lightweight S3 objects (named lists). Population references may be
# string names or zero-based integer IDs; they are resolved against the
# Demography object when the simulation is set up.

.make_event <- function(cls, super, fields) {
  structure(fields, class = c(cls, super, "DemographicEvent", "msprimeR_object"))
}

.check_event_time <- function(time) {
  if (!.is_number(time)) stop("Event time must be a single number", call. = FALSE)
  as.double(time)
}

#' Demographic events
#'
#' Constructors for the demographic events of the msprime 0.x API, retained for
#' use with [simulate()] and `Demography$from_old_style()`. In new code, use
#' the corresponding `Demography` methods (e.g. `add_population_split()`,
#' `add_mass_migration()`, `add_population_parameters_change()`).
#'
#' Population references are zero-based integer IDs or population names.
#'
#' @param time The time (in generations ago) at which the event occurs.
#' @param initial_size The new population size.
#' @param growth_rate The new exponential growth rate.
#' @param population The population affected (`NULL` or -1 for all).
#' @param population_id Deprecated alias for `population`.
#' @param rate The new migration rate.
#' @param source,dest Source and destination populations (backwards in time).
#' @param matrix_index Deprecated: a `c(source, dest)` pair.
#' @param proportion The probability that each lineage moves (mass migration)
#'   or coalesces (simple bottleneck).
#' @param destination Deprecated alias for `dest`.
#' @param strength Equivalent duration of the coalescent process.
#' @return A `DemographicEvent` object.
#' @name demographic_events
NULL

#' @rdname demographic_events
#' @export
PopulationParametersChange <- function(time, initial_size = NULL, growth_rate = NULL,
                                       population = NULL, population_id = NULL) {
  if (!is.null(population_id) && !is.null(population)) {
    stop("population_id and population are aliases; cannot supply both.", call. = FALSE)
  }
  if (!is.null(population_id)) population <- population_id
  if (is.null(growth_rate) && is.null(initial_size)) {
    stop("Must specify one or more of growth_rate and initial_size", call. = FALSE)
  }
  if (!is.null(initial_size) && initial_size < 0) {
    stop("Cannot have a population size < 0", call. = FALSE)
  }
  .make_event("PopulationParametersChange", "ParameterChangeEvent", list(
    time = .check_event_time(time), initial_size = initial_size,
    growth_rate = growth_rate, population = if (is.null(population)) -1 else population))
}

#' @rdname demographic_events
#' @export
MigrationRateChange <- function(time, rate, source = -1, dest = -1, matrix_index = NULL) {
  if (!is.null(matrix_index)) {
    source <- matrix_index[[1]]
    dest <- matrix_index[[2]]
  }
  .make_event("MigrationRateChange", "ParameterChangeEvent", list(
    time = .check_event_time(time), rate = rate, source = source, dest = dest))
}

#' @rdname demographic_events
#' @export
MassMigration <- function(time, source, dest = NULL, proportion = 1.0, destination = NULL) {
  if (!is.null(dest) && !is.null(destination)) {
    stop("dest and destination are aliases; cannot supply both", call. = FALSE)
  }
  if (!is.null(destination)) dest <- destination
  .make_event("MassMigration", "LineageMovementEvent", list(
    time = .check_event_time(time), source = source, dest = dest, proportion = proportion))
}

#' @rdname demographic_events
#' @export
CensusEvent <- function(time) {
  .make_event("CensusEvent", character(0), list(time = .check_event_time(time)))
}

#' @rdname demographic_events
#' @export
SimpleBottleneck <- function(time, population, proportion = 1.0) {
  .make_event("SimpleBottleneck", "StateChangeEvent", list(
    time = .check_event_time(time), population = population, proportion = proportion))
}

#' @rdname demographic_events
#' @export
InstantaneousBottleneck <- function(time, population, strength = 1.0) {
  .make_event("InstantaneousBottleneck", "StateChangeEvent", list(
    time = .check_event_time(time), population = population, strength = strength))
}

.SymmetricMigrationRateChange <- function(time, populations, rate) {
  .make_event("SymmetricMigrationRateChange", "ParameterChangeEvent", list(
    time = .check_event_time(time), populations = as.list(populations), rate = rate))
}

.ActivatePopulationEvent <- function(time, population) {
  .make_event("ActivatePopulationEvent", character(0), list(
    time = .check_event_time(time), population = population))
}

.PopulationSplit <- function(time, derived, ancestral) {
  .make_event("PopulationSplit", "LineageMovementEvent", list(
    time = .check_event_time(time), derived = as.list(derived), ancestral = ancestral))
}

.Admixture <- function(time, derived, ancestral, proportions) {
  .make_event("Admixture", "LineageMovementEvent", list(
    time = .check_event_time(time), derived = derived, ancestral = as.list(ancestral),
    proportions = as.double(unlist(proportions))))
}

.event_type_str <- c(
  PopulationParametersChange = "Population parameter change",
  MigrationRateChange = "Migration rate change",
  SymmetricMigrationRateChange = "Symmetric migration rate change",
  MassMigration = "Mass Migration",
  ActivatePopulationEvent = "Activate population event",
  PopulationSplit = "Population Split",
  Admixture = "Admixture",
  SimpleBottleneck = "Simple Bottleneck",
  InstantaneousBottleneck = "Instantaneous Bottleneck",
  CensusEvent = "Census")

# Resolve a population reference into a zero-based integer ID (-1 means all).
.convert_id <- function(ref, demography) {
  if (is.null(ref)) return(-1L)
  if (is.numeric(ref) && length(ref) == 1 && !is.na(ref) && ref == -1) return(-1L)
  if (is.null(demography)) {
    if (!.isinteger(ref)) {
      stop("Working with demographic events not associated with a Demography object ",
           "is a legacy-only operation. Population references must be integer IDs",
           call. = FALSE)
    }
    return(as.integer(ref))
  }
  as.integer(demography$get(ref)$id)
}

# Low-level representation passed to the C simulator.
.ll_event <- function(event, demography) {
  cls <- class(event)[1]
  conv <- function(ref) .convert_id(ref, demography)
  switch(cls,
    PopulationParametersChange = {
      ret <- list(type = "population_parameters_change", time = event$time,
                  population = conv(event$population))
      if (!is.null(event$growth_rate)) ret$growth_rate <- event$growth_rate
      if (!is.null(event$initial_size)) ret$initial_size <- event$initial_size
      ret
    },
    MigrationRateChange = list(type = "migration_rate_change", time = event$time,
                               migration_rate = event$rate, source = conv(event$source),
                               dest = conv(event$dest)),
    SymmetricMigrationRateChange = list(
      type = "symmetric_migration_rate_change", time = event$time,
      populations = vapply(event$populations, conv, integer(1)), rate = event$rate),
    MassMigration = list(type = "mass_migration", time = event$time,
                         source = conv(event$source), dest = conv(event$dest),
                         proportion = event$proportion),
    ActivatePopulationEvent = list(type = "activate_population_event", time = event$time,
                                   population = conv(event$population)),
    PopulationSplit = list(type = "population_split", time = event$time,
                           derived = vapply(event$derived, conv, integer(1)),
                           ancestral = conv(event$ancestral)),
    Admixture = list(type = "admixture", time = event$time, derived = conv(event$derived),
                     ancestral = vapply(event$ancestral, conv, integer(1)),
                     proportions = event$proportions),
    SimpleBottleneck = list(type = "simple_bottleneck", time = event$time,
                            population = conv(event$population),
                            proportion = event$proportion),
    InstantaneousBottleneck = list(type = "instantaneous_bottleneck", time = event$time,
                                   population = conv(event$population),
                                   strength = event$strength),
    CensusEvent = list(type = "census_event", time = event$time),
    .stopf("Unknown demographic event class '%s'", cls))
}

# Conditional ("sequential") proportions <-> absolute proportions.
.sequential_to_proportions <- function(S) {
  P <- numeric(0)
  for (j in seq_along(S)) P <- c(P, S[j] * (1 - .psum(P)))
  P
}

.proportions_to_sequential <- function(P) {
  C <- numeric(length(P))
  for (j in seq_along(P)) {
    s <- .psum(P[seq_len(j - 1)])
    if (s < 1) C[j] <- P[j] / (1 - s)
  }
  C
}

# Lineage movements (source, dest, proportion) as zero-based IDs.
.as_lineage_movements <- function(event, demography) {
  conv <- function(ref) .convert_id(ref, demography)
  cls <- class(event)[1]
  if (cls == "MassMigration") {
    return(list(list(source = conv(event$source), dest = conv(event$dest),
                     proportion = event$proportion)))
  }
  if (cls == "PopulationSplit") {
    anc <- conv(event$ancestral)
    return(lapply(event$derived, function(pop) {
      list(source = conv(pop), dest = anc, proportion = 1)
    }))
  }
  if (cls == "Admixture") {
    derived <- conv(event$derived)
    S <- .proportions_to_sequential(event$proportions)
    return(lapply(seq_along(event$ancestral), function(j) {
      list(source = derived, dest = conv(event$ancestral[[j]]), proportion = S[j])
    }))
  }
  stop("Not a lineage movement event", call. = FALSE)
}

.list_str <- function(a, fmt = NULL) {
  a <- unlist(a)
  items <- if (is.null(fmt)) format(a, trim = TRUE) else sprintf(fmt, a)
  paste0("[", paste(items, collapse = ", "), "]")
}

.g <- function(x, digits) formatC(x, digits = digits, format = "g")

.event_parameters <- function(event) {
  switch(class(event)[1],
    PopulationParametersChange = {
      s <- sprintf("population=%s", event$population)
      if (!is.null(event$initial_size)) s <- paste0(s, ", initial_size=", event$initial_size)
      if (!is.null(event$growth_rate)) s <- paste0(s, ", growth_rate=", event$growth_rate)
      s
    },
    MigrationRateChange = sprintf("source=%s, dest=%s, rate=%s", event$source,
                                  event$dest, event$rate),
    SymmetricMigrationRateChange = sprintf("populations=%s, rate=%s",
                                           .list_str(event$populations), event$rate),
    MassMigration = sprintf("source=%s, dest=%s, proportion=%s", event$source,
                            event$dest, .g(event$proportion, 3)),
    ActivatePopulationEvent = sprintf("population=%s", event$population),
    PopulationSplit = sprintf("derived=%s, ancestral=%s", .list_str(event$derived),
                              event$ancestral),
    Admixture = sprintf("derived=%s ancestral=%s proportions=%s", event$derived,
                        .list_str(event$ancestral), .list_str(event$proportions, "%.2f")),
    SimpleBottleneck = sprintf("population=%s, proportion=%s", event$population,
                               event$proportion),
    InstantaneousBottleneck = sprintf("population=%s, strength=%s", event$population,
                                      event$strength),
    CensusEvent = "")
}

.event_effect <- function(event) {
  switch(class(event)[1],
    PopulationParametersChange = {
      s <- ""
      if (!is.null(event$initial_size)) {
        s <- paste0(s, "initial_size \u2192 ", .g(event$initial_size, 2), " ")
        if (!is.null(event$growth_rate)) s <- paste0(s, "and ")
      }
      if (!is.null(event$growth_rate)) {
        s <- paste0(s, "growth_rate \u2192 ", .g(event$growth_rate, 3), " ")
      }
      s <- paste0(s, "for")
      if (identical(as.double(event$population), -1)) {
        paste0(s, " all populations")
      } else {
        paste0(s, " population ", event$population)
      }
    },
    MigrationRateChange = {
      ret <- "Backwards-time migration rate "
      if (identical(as.double(event$source), -1) && identical(as.double(event$dest), -1)) {
        ret <- paste0(ret, "for all populations ")
      } else {
        ret <- paste0(ret, sprintf("from %s to %s ", event$source, event$dest))
      }
      paste0(ret, "\u2192 ", event$rate)
    },
    SymmetricMigrationRateChange = {
      s <- "Sets the symmetric migration rate between "
      if (length(event$populations) == 2) {
        s <- paste0(s, sprintf("%s and %s ", event$populations[[1]], event$populations[[2]]))
      } else {
        s <- paste0(s, sprintf("all pairs of populations in %s ",
                               .list_str(event$populations)))
      }
      paste0(s, sprintf("to %s per generation", event$rate))
    },
    MassMigration = {
      ret <- if (event$proportion == 1) {
        sprintf("All lineages currently in population %s move to %s ", event$source,
                event$dest)
      } else {
        sprintf("Lineages currently in population %s move to %s with probability %s ",
                event$source, event$dest, .g(event$proportion, 3))
      }
      paste0(ret, sprintf("(equivalent to individuals migrating from %s to %s forwards in time)",
                          event$dest, event$source))
    },
    ActivatePopulationEvent = sprintf("Activates population %s", event$population),
    PopulationSplit = {
      d <- event$derived
      s <- "Moves all lineages from "
      if (length(d) == 1) {
        s <- paste0(s, sprintf("the '%s' derived population ", d[[1]]))
      } else if (length(d) == 2) {
        s <- paste0(s, sprintf("derived populations '%s' and '%s' ", d[[1]], d[[2]]))
      } else {
        s <- paste0(s, sprintf("derived populations %s ", .list_str(d)))
      }
      s <- paste0(s, sprintf("to the ancestral '%s' population. Also set ", event$ancestral))
      s <- paste0(s, if (length(d) == 1) sprintf("'%s' ", d[[1]]) else "the derived populations ")
      paste0(s, "to inactive, and all migration rates to and from the derived population",
             if (length(d) > 1) "s" else "", " to zero.")
    },
    Admixture = {
      move_to <- paste(vapply(seq_along(event$ancestral), function(j) {
        sprintf("'%s' with proba %s", event$ancestral[[j]], .g(event$proportions[j], 3))
      }, character(1)), collapse = "; ")
      sprintf(paste0("Moves all lineages from admixed population '%s' to ancestral ",
                     "population%s. Lineages move to %s. Set '%s' to inactive, and all ",
                     "migration rates to and from '%s' to zero."),
              event$derived, if (length(event$ancestral) > 1) "s" else "", move_to,
              event$derived, event$derived)
    },
    SimpleBottleneck = sprintf("Lineages in population %s coalesce with probability %s",
                               event$population, event$proportion),
    InstantaneousBottleneck = sprintf("Equivalent to %s generations of the coalescent",
                                      event$strength),
    CensusEvent = "Insert census nodes to record the location of all lineages")
}

#' @export
print.DemographicEvent <- function(x, ...) {
  cat(sprintf("%s(time=%s, %s)\n", class(x)[1], format(x$time), .event_parameters(x)))
  invisible(x)
}

#' @export
format.DemographicEvent <- function(x, ...) {
  sprintf("%s(time=%s, %s)", class(x)[1], format(x$time), .event_parameters(x))
}

.events_equal <- function(e1, e2) {
  identical(class(e1), class(e2)) && isTRUE(all.equal(unclass(e1), unclass(e2),
                                                       tolerance = 0))
}
