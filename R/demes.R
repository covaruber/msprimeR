# Demes (https://popsim-consortium.github.io/demes-spec-docs/) support:
# resolution of the demes data model, and conversion to/from Demography.

.DEMES_TOL <- list(rel = 1e-9, abs = 1e-12)

.as_time <- function(x, what) {
  if (is.character(x) && length(x) == 1 && x %in% c("Infinity", "inf", ".inf")) return(Inf)
  x <- .as_num(x)
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) .stopf("%s must be a number", what)
  as.double(x)
}

# The 'yaml' package reads numbers such as 1e-4 (valid YAML 1.2, but not 1.1)
# as strings.
.as_num <- function(x) {
  if (is.list(x)) x <- unlist(x)
  if (is.character(x)) {
    v <- suppressWarnings(as.double(x))
    if (!anyNA(v)) return(v)
  }
  x
}

.check_fields <- function(data, allowed, scope) {
  bad <- setdiff(names(data), allowed)
  if (length(bad) > 0) {
    .stopf("%s: unexpected field: '%s'. Allowed fields are: %s", scope, bad[1],
           paste(allowed, collapse = ", "))
  }
}

.insert_defaults <- function(data, defaults) {
  for (k in names(defaults)) if (is.null(data[[k]])) data[[k]] <- defaults[[k]]
  data
}

.demes_input <- function(graph) {
  if (inherits(graph, "demes_graph")) return(graph)
  if (is.character(graph) && length(graph) == 1) {
    return(resolve_demes(.read_yaml(graph)))
  }
  if (is.list(graph)) return(resolve_demes(graph))
  stop("graph must be a demes graph, a list, a YAML string or a file path", call. = FALSE)
}

# Read YAML with YAML 1.2 booleans (only true/false), as used by demes: the
# 'yaml' package otherwise follows YAML 1.1, where e.g. a deme named "Y" or
# "no" would be read as a logical.
.read_yaml <- function(x) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The 'yaml' package is required to read demes models", call. = FALSE)
  }
  handlers <- list(
    "bool#yes" = function(v) if (v %in% c("true", "True", "TRUE")) TRUE else v,
    "bool#no" = function(v) if (v %in% c("false", "False", "FALSE")) FALSE else v)
  if (file.exists(x)) yaml::read_yaml(x, handlers = handlers) else
    yaml::yaml.load(x, handlers = handlers)
}

#' Resolve a demes model
#'
#' Resolves and validates a demes model given in the "human data model" (for
#' example as read from YAML by `yaml::read_yaml()`), filling in all default
#' values, as `demes.Graph.fromdict()` does in Python.
#'
#' @param data A named list following the demes specification.
#' @return A fully resolved demes graph (class `demes_graph`).
#' @export
resolve_demes <- function(data) {
  if (inherits(data, "demes_graph")) return(data)
  if (!is.list(data)) stop("data is not a dictionary", call. = FALSE)
  .check_fields(data, c("description", "time_units", "generation_time", "defaults", "doi",
                        "metadata", "demes", "migrations", "pulses"), "toplevel")
  defaults <- data$defaults %||% list()
  .check_fields(defaults, c("deme", "migration", "pulse", "epoch"), "defaults")
  deme_defaults <- defaults$deme %||% list()
  migration_defaults <- defaults$migration %||% list()
  pulse_defaults <- defaults$pulse %||% list()
  global_epoch_defaults <- defaults$epoch %||% list()
  if (is.null(data$time_units)) stop("toplevel: required field 'time_units' not found",
                                     call. = FALSE)
  graph <- list(description = data$description %||% "", time_units = data$time_units,
                generation_time = .as_num(data$generation_time),
                doi = as.list(data$doi %||% list()),
                metadata = data$metadata %||% list(), demes = list(), migrations = list(),
                pulses = list())
  if (graph$time_units != "generations" && is.null(graph$generation_time)) {
    stop("if time_units!=\"generations\", generation_time must be specified", call. = FALSE)
  }
  if (is.null(graph$generation_time)) graph$generation_time <- 1
  if (graph$time_units == "generations" && graph$generation_time != 1) {
    stop("time_units==\"generations\", but generation_time!=1", call. = FALSE)
  }
  demes_list <- data$demes
  if (length(demes_list) == 0) stop("toplevel: 'demes' must be a non-empty list",
                                    call. = FALSE)
  deme_index <- function(name) {
    for (j in seq_along(graph$demes)) if (graph$demes[[j]]$name == name) return(j)
    NA_integer_
  }
  for (i in seq_along(demes_list)) {
    dd <- demes_list[[i]]
    if (is.null(dd$name)) .stopf("demes[%d]: required field 'name' not found", i - 1)
    name <- dd$name
    scope <- sprintf("demes[%d] %s", i - 1, name)
    .check_fields(dd, c("name", "description", "start_time", "ancestors", "proportions",
                        "defaults", "epochs"), scope)
    if (!.is_identifier(name)) .stopf("Invalid deme name '%s'", name)
    if (!is.na(deme_index(name))) .stopf("%s: field 'name' must be unique", scope)
    dd <- .insert_defaults(dd, deme_defaults)
    ancestors <- as.character(unlist(dd$ancestors %||% list()))
    for (a in ancestors) {
      if (is.na(deme_index(a))) {
        .stopf(paste0("%s: ancestor deme '%s' not found. Note: ancestor demes must be ",
                      "specified before their children."), scope, a)
      }
    }
    proportions <- as.double(.as_num(dd$proportions))
    if (is.null(dd$proportions)) proportions <- if (length(ancestors) == 1) 1 else numeric(0)
    if (is.null(dd$start_time)) {
      if (length(ancestors) > 1) {
        .stopf(paste0("%s: field 'start_time' not found, but is required for demes with ",
                      "multiple ancestors"), scope)
      }
      start_time <- if (length(ancestors) == 1) {
        anc <- graph$demes[[deme_index(ancestors)]]
        anc$epochs[[length(anc$epochs)]]$end_time
      } else {
        Inf
      }
    } else {
      start_time <- .as_time(dd$start_time, paste(scope, "start_time"))
    }
    if (length(ancestors) == 0 && !is.infinite(start_time)) {
      .stopf(paste0("%s: field 'ancestors' not found, but is required for demes with a ",
                    "finite 'start_time'"), scope)
    }
    for (a in ancestors) {
      anc <- graph$demes[[deme_index(a)]]
      anc_end <- anc$epochs[[length(anc$epochs)]]$end_time
      if (!(anc$start_time > start_time && start_time >= anc_end)) {
        .stopf(paste0("%s: start_time=%s is outside the interval of existence for ancestor ",
                      "'%s'"), scope, start_time, a)
      }
    }
    if (anyDuplicated(ancestors)) .stopf("deme %s: duplicate ancestors", name)
    if (name %in% ancestors) .stopf("deme %s: deme cannot be its own ancestor", name)
    if (length(ancestors) != length(proportions)) {
      .stopf("deme %s: ancestors and proportions have different lengths", name)
    }
    if (length(proportions) > 0 && !.isclose(sum(proportions), 1)) {
      .stopf("deme %s: ancestry proportions must sum to 1.0", name)
    }
    local_defaults <- dd$defaults %||% list()
    .check_fields(local_defaults, "epoch", paste(scope, "defaults"))
    epoch_defaults <- global_epoch_defaults
    for (k in names(local_defaults$epoch)) epoch_defaults[[k]] <- local_defaults$epoch[[k]]
    if (length(epoch_defaults) == 0 && is.null(dd$epochs)) {
      .stopf("%s: required field 'epochs' not found", scope)
    }
    epochs_in <- dd$epochs %||% list(list())
    if (length(epochs_in) == 0) .stopf("%s: 'epochs' must be a non-empty list", scope)
    epochs <- list()
    for (j in seq_along(epochs_in)) {
      ed <- epochs_in[[j]] %||% list()
      escope <- sprintf("%s: epochs[%d]", scope, j - 1)
      .check_fields(ed, c("end_time", "start_size", "end_size", "size_function",
                          "cloning_rate", "selfing_rate"), escope)
      ed <- .insert_defaults(ed, epoch_defaults)
      if (is.null(ed$end_time)) {
        if (j == length(epochs_in)) ed$end_time <- 0 else
          .stopf("%s: required field 'end_time' not found", escope)
      }
      for (f in c("start_size", "end_size", "selfing_rate", "cloning_rate")) {
        if (!is.null(ed[[f]])) ed[[f]] <- .as_num(ed[[f]])
      }
      end_time <- .as_time(ed$end_time, paste(escope, "end_time"))
      if (length(epochs) == 0) {
        e_start <- start_time
        if (is.null(ed$start_size) && is.null(ed$end_size)) {
          .stopf("deme %s: first epoch must have start_size or end_size", name)
        }
        start_size <- ed$start_size %||% ed$end_size
        end_size <- ed$end_size %||% start_size
      } else {
        prev <- epochs[[length(epochs)]]
        e_start <- prev$end_time
        start_size <- ed$start_size %||% prev$end_size
        end_size <- ed$end_size %||% start_size
      }
      size_function <- ed$size_function %||%
        if (start_size == end_size) "constant" else "exponential"
      epoch <- list(start_time = e_start, end_time = end_time,
                    start_size = as.double(start_size), end_size = as.double(end_size),
                    size_function = size_function,
                    selfing_rate = as.double(ed$selfing_rate %||% 0),
                    cloning_rate = as.double(ed$cloning_rate %||% 0))
      if (!(epoch$start_time > epoch$end_time)) {
        .stopf("deme %s: epoch[%d]: invalid epoch: must have start_time > end_time", name,
               j - 1)
      }
      if (epoch$end_time < 0 || is.infinite(epoch$end_time)) {
        .stopf("deme %s: epoch[%d]: invalid end_time", name, j - 1)
      }
      if (!(epoch$start_size > 0 && epoch$end_size > 0) || !is.finite(epoch$start_size) ||
          !is.finite(epoch$end_size)) {
        .stopf("deme %s: epoch[%d]: sizes must be positive and finite", name, j - 1)
      }
      if (!size_function %in% c("constant", "exponential", "linear")) {
        .stopf("deme %s: epoch[%d]: invalid size_function", name, j - 1)
      }
      if (is.infinite(epoch$start_time) && epoch$start_size != epoch$end_size) {
        .stopf("deme %s: epoch[%d]: if start time is inf, must be a constant size epoch",
               name, j - 1)
      }
      if (size_function == "constant" && epoch$start_size != epoch$end_size) {
        .stopf("deme %s: epoch[%d]: start_size != end_size, but size_function is constant",
               name, j - 1)
      }
      epochs[[length(epochs) + 1]] <- epoch
    }
    graph$demes[[length(graph$demes) + 1]] <- list(
      name = name, description = dd$description %||% "", start_time = start_time,
      ancestors = ancestors, proportions = proportions, epochs = epochs)
  }
  deme_interval <- function(name) {
    d <- graph$demes[[deme_index(name)]]
    c(start = d$start_time, end = d$epochs[[length(d$epochs)]]$end_time)
  }
  intersection <- function(a, b, time = NULL) {
    ia <- deme_interval(a)
    ib <- deme_interval(b)
    lo <- max(ia[["end"]], ib[["end"]])
    hi <- min(ia[["start"]], ib[["start"]])
    if (!is.null(time) && !(lo <= time && time <= hi)) {
      .stopf("time %s not in interval [%s, %s], as defined by the time-intersection of %s and %s",
             time, lo, hi, a, b)
    }
    c(lo = lo, hi = hi)
  }
  add_asymmetric <- function(source, dest, rate, start_time, end_time) {
    for (d in c(source, dest)) if (is.na(deme_index(d))) .stopf("%s not in graph", d)
    if (source == dest) stop("source and dest cannot be the same deme", call. = FALSE)
    iv <- intersection(source, dest, start_time)
    start_time <- start_time %||% iv[["hi"]]
    end_time <- end_time %||% iv[["lo"]]
    intersection(source, dest, end_time)
    if (!(start_time > end_time)) stop("must have start_time > end_time", call. = FALSE)
    if (!(rate >= 0 && rate <= 1)) stop("must have 0 <= rate <= 1", call. = FALSE)
    graph$migrations[[length(graph$migrations) + 1]] <<- list(
      source = source, dest = dest, start_time = start_time, end_time = end_time,
      rate = as.double(rate))
  }
  for (i in seq_along(data$migrations)) {
    md <- .insert_defaults(data$migrations[[i]], migration_defaults)
    .check_fields(md, c("demes", "source", "dest", "start_time", "end_time", "rate"),
                  sprintf("migration[%d]", i - 1))
    if (is.null(md$rate)) .stopf("migration[%d]: required field 'rate' not found", i - 1)
    md$rate <- .as_num(md$rate)
    st <- if (is.null(md$start_time)) NULL else .as_time(md$start_time, "start_time")
    et <- if (is.null(md$end_time)) NULL else .as_time(md$end_time, "end_time")
    if (!is.null(md$demes) && is.null(md$source) && is.null(md$dest)) {
      demes <- as.character(unlist(md$demes))
      if (length(demes) < 2) stop("must specify a list of two or more deme names",
                                  call. = FALSE)
      for (a in demes) for (b in demes) if (a != b) add_asymmetric(a, b, md$rate, st, et)
    } else if (is.null(md$demes) && !is.null(md$source) && !is.null(md$dest)) {
      add_asymmetric(md$source, md$dest, md$rate, st, et)
    } else {
      .stopf(paste0("migration[%d]: must be symmetric (specify 'demes' list) *or* ",
                    "asymmetric (specify both 'source' and 'dest')."), i - 1)
    }
  }
  for (i in seq_along(data$pulses)) {
    pd <- .insert_defaults(data$pulses[[i]], pulse_defaults)
    .check_fields(pd, c("sources", "dest", "time", "proportions"), sprintf("pulse[%d]", i - 1))
    for (f in c("sources", "dest", "time", "proportions")) {
      if (is.null(pd[[f]])) .stopf("pulse[%d]: required field '%s' not found", i - 1, f)
    }
    sources <- as.character(unlist(pd$sources))
    proportions <- as.double(.as_num(pd$proportions))
    time <- .as_time(pd$time, "time")
    for (d in c(sources, pd$dest)) if (is.na(deme_index(d))) .stopf("%s not in graph", d)
    for (s in sources) intersection(s, pd$dest, time)
    if (time == deme_interval(pd$dest)[["end"]]) {
      .stopf("invalid pulse at time=%s, which is dest=%s's end_time", time, pd$dest)
    }
    for (s in sources) {
      if (time == deme_interval(s)[["start"]]) {
        .stopf("invalid pulse at time=%s, which is source=%s's start_time", time, s)
      }
    }
    if (length(sources) != length(proportions)) {
      stop("sources and proportions must have the same length", call. = FALSE)
    }
    if (sum(proportions) > 1 + 1e-12) stop("proportions must sum to less than one",
                                           call. = FALSE)
    graph$pulses[[length(graph$pulses) + 1]] <- list(sources = sources, dest = pd$dest,
                                                     time = time, proportions = proportions)
  }
  if (length(graph$pulses) > 1) {
    times <- vapply(graph$pulses, function(p) p$time, numeric(1))
    graph$pulses <- graph$pulses[order(-times, method = "radix")]
  }
  structure(graph, class = "demes_graph")
}

#' @export
print.demes_graph <- function(x, ...) {
  cat(sprintf("<demes_graph> time_units: %s, %d demes, %d migrations, %d pulses\n",
              x$time_units, length(x$demes), length(x$migrations), length(x$pulses)))
  for (d in x$demes) {
    cat(sprintf("  %s: start_time=%s, ancestors=[%s], %d epoch(s)\n", d$name,
                format(d$start_time), paste(d$ancestors, collapse = ", "), length(d$epochs)))
  }
  invisible(x)
}

.demes_in_generations <- function(graph) {
  g <- graph$generation_time
  for (i in seq_along(graph$demes)) {
    graph$demes[[i]]$start_time <- graph$demes[[i]]$start_time / g
    for (j in seq_along(graph$demes[[i]]$epochs)) {
      graph$demes[[i]]$epochs[[j]]$start_time <- graph$demes[[i]]$epochs[[j]]$start_time / g
      graph$demes[[i]]$epochs[[j]]$end_time <- graph$demes[[i]]$epochs[[j]]$end_time / g
    }
  }
  for (i in seq_along(graph$migrations)) {
    graph$migrations[[i]]$start_time <- graph$migrations[[i]]$start_time / g
    graph$migrations[[i]]$end_time <- graph$migrations[[i]]$end_time / g
  }
  for (i in seq_along(graph$pulses)) graph$pulses[[i]]$time <- graph$pulses[[i]]$time / g
  graph$time_units <- "generations"
  graph$generation_time <- 1
  graph
}

.deme_end_time <- function(d) d$epochs[[length(d$epochs)]]$end_time

.demes_discrete_events <- function(graph) {
  by_name <- stats::setNames(graph$demes, vapply(graph$demes, function(d) d$name, ""))
  ev <- list(pulses = graph$pulses, splits = list(), branches = list(), mergers = list(),
             admixtures = list())
  splits_to_add <- list()
  for (d in graph$demes) {
    p <- d$ancestors
    if (length(p) == 0) next
    if (length(p) == 1) {
      if (d$start_time == .deme_end_time(by_name[[p]])) {
        splits_to_add[[p]] <- c(splits_to_add[[p]], d$name)
      } else {
        ev$branches[[length(ev$branches) + 1]] <- list(parent = p, child = d$name,
                                                       time = d$start_time)
      }
    } else {
      aligned <- all(vapply(p, function(a) d$start_time == .deme_end_time(by_name[[a]]),
                            logical(1)))
      item <- list(parents = p, proportions = d$proportions, child = d$name,
                   time = d$start_time)
      if (aligned) ev$mergers[[length(ev$mergers) + 1]] <- item else
        ev$admixtures[[length(ev$admixtures) + 1]] <- item
    }
  }
  for (parent in names(splits_to_add)) {
    ev$splits[[length(ev$splits) + 1]] <- list(parent = parent,
                                               children = splits_to_add[[parent]],
                                               time = .deme_end_time(by_name[[parent]]))
  }
  ev
}

.demography_from_demes <- function(graph) {
  graph <- .demes_input(graph)
  for (d in graph$demes) {
    for (e in d$epochs) {
      if (e$selfing_rate != 0) .stopf("deme %s: non-zero selfing_rate not supported", d$name)
      if (e$cloning_rate != 0) .stopf("deme %s: non-zero cloning_rate not supported", d$name)
    }
  }
  get_growth_rate <- function(epoch) {
    if (!epoch$size_function %in% c("constant", "exponential")) {
      stop("msprime only supports constant or exponentially changing population sizes",
           call. = FALSE)
    }
    if (epoch$end_size != epoch$start_size) {
      -log(epoch$start_size / epoch$end_size) / (epoch$start_time - epoch$end_time)
    } else {
      0
    }
  }
  graph <- .demes_in_generations(graph)
  events <- .demes_discrete_events(graph)
  demography <- Demography$new()
  for (d in graph$demes) {
    initial_size <- 0
    growth_rate <- 0
    last_epoch <- d$epochs[[length(d$epochs)]]
    initially_active <- FALSE
    to_activate <- FALSE
    if (last_epoch$end_time == 0) {
      initial_size <- last_epoch$end_size
      growth_rate <- get_growth_rate(last_epoch)
      initially_active <- TRUE
    } else {
      to_activate <- TRUE
      for (s in events$splits) if (d$name == s$parent) to_activate <- FALSE
      for (m in events$mergers) if (d$name %in% m$parents) to_activate <- FALSE
    }
    demography$add_population(name = d$name, description = d$description,
                              growth_rate = growth_rate, initial_size = initial_size,
                              default_sampling_time = .deme_end_time(d),
                              initially_active = initially_active)
    if (to_activate) {
      demography$add_activate_population_event(time = .deme_end_time(d), population = d$name)
    }
    for (epoch in rev(d$epochs)) {
      new_initial_size <- if (initial_size != epoch$end_size) epoch$end_size else NULL
      alpha <- get_growth_rate(epoch)
      new_growth_rate <- if (growth_rate != alpha) alpha else NULL
      if (!is.null(new_growth_rate) || !is.null(new_initial_size)) {
        demography$add_population_parameters_change(
          population = d$name, time = epoch$end_time, initial_size = new_initial_size,
          growth_rate = new_growth_rate)
        # msprime resets both tracked values here (even the unchanged one),
        # which adds redundant events for later epochs; track the true values.
        if (!is.null(new_initial_size)) initial_size <- new_initial_size
        if (!is.null(new_growth_rate)) growth_rate <- new_growth_rate
      }
    }
  }
  for (s in events$splits) {
    demography$add_population_split(time = s$time, ancestral = s$parent, derived = s$children)
  }
  for (b in events$branches) {
    demography$add_population_split(time = b$time, ancestral = b$parent, derived = b$child)
  }
  for (pulse in rev(events$pulses)) {
    seq_props <- .proportions_to_sequential(pulse$proportions)
    for (j in seq_along(pulse$sources)) {
      demography$add_mass_migration(time = pulse$time, source = pulse$dest,
                                    dest = pulse$sources[j], proportion = seq_props[j])
    }
  }
  for (m in c(events$mergers, events$admixtures)) {
    demography$add_admixture(time = m$time, derived = m$child, ancestral = m$parents,
                             proportions = m$proportions)
  }
  for (mig in graph$migrations) {
    if (!is.infinite(mig$start_time)) {
      demography$add_migration_rate_change(time = mig$start_time, source = mig$dest,
                                           dest = mig$source, rate = 0)
    }
  }
  for (mig in graph$migrations) {
    if (mig$end_time == 0) {
      demography$set_migration_rate(source = mig$dest, dest = mig$source, rate = mig$rate)
    } else {
      demography$add_migration_rate_change(time = mig$end_time, source = mig$dest,
                                           dest = mig$source, rate = mig$rate)
    }
  }
  demography$sort_events()
  demography
}

.demography_to_demes <- function(demography) {
  dbg <- demography$debug()
  resolved <- dbg$demography
  pop_names <- names(resolved)
  demes <- list()
  for (pop in resolved$populations) {
    active <- Filter(function(e) e$populations[[pop$id + 1]]$active, dbg$epochs)
    start_time <- max(vapply(active, function(e) e$end_time, numeric(1)))
    end_time <- min(vapply(active, function(e) e$start_time, numeric(1)))
    demes[[pop$name]] <- list(
      name = pop$name,
      description = if (nzchar(pop$description %||% "")) pop$description else NULL,
      start_time = start_time,
      epochs = list(list(end_size = pop$initial_size, growth_rate = pop$growth_rate,
                         end_time = end_time)))
  }
  epoch_resolve <- function(name, time) {
    deme <- demes[[name]]
    epoch <- deme$epochs[[1]]
    if (!(deme$start_time > time && time >= epoch$end_time)) {
      .stopf("time %s outside %s's existence interval (start_time=%s, end_time=%s]",
             time, name, deme$start_time, epoch$end_time)
    }
    if (time > epoch$end_time) {
      new_epoch <- epoch
      growth_rate <- epoch$growth_rate %||% 0
      epoch$growth_rate <- NULL
      dt <- time - epoch$end_time
      size_at_t <- epoch$end_size * exp(-growth_rate * dt)
      epoch$start_size <- size_at_t
      new_epoch$end_size <- size_at_t
      new_epoch$end_time <- time
      deme$epochs[[1]] <- epoch
      deme$epochs <- c(list(new_epoch), deme$epochs)
    }
    demes[[name]] <<- deme
  }
  pulses <- list()
  events <- resolved$events
  times <- vapply(events, function(e) e$time, numeric(1))
  groups <- split(seq_along(events), factor(times, levels = unique(times)))
  for (g in groups) {
    time <- events[[g[1]]]$time
    lineage_movements <- list()
    for (event in events[g]) {
      if (inherits(event, "LineageMovementEvent")) {
        lineage_movements <- c(lineage_movements, .as_lineage_movements(event, resolved))
      } else if (inherits(event, "PopulationParametersChange")) {
        pids <- if (identical(as.double(event$population), -1)) {
          seq_len(resolved$num_populations) - 1
        } else {
          resolved$get(event$population)$id
        }
        for (pid in pids) {
          name <- pop_names[pid + 1]
          deme <- demes[[name]]
          if (length(deme$epochs) == 1 && deme$epochs[[1]]$end_size == 0 &&
              !is.null(event$initial_size) && event$initial_size > 0) {
            demes[[name]]$epochs[[1]]$end_time <- time
          }
          epoch_resolve(name, time)
          if (!is.null(event$growth_rate)) demes[[name]]$epochs[[1]]$growth_rate <- event$growth_rate
          if (!is.null(event$initial_size)) demes[[name]]$epochs[[1]]$end_size <- event$initial_size
        }
      } else if (inherits(event, "StateChangeEvent")) {
        .stopf("Cannot convert %s to demes model.", format(event))
      }
    }
    by_source <- list()
    for (lm in lineage_movements) {
      key <- pop_names[lm$source + 1]
      by_source[[key]] <- c(by_source[[key]], list(lm))
    }
    pulse_keys <- character(0)
    for (source in names(by_source)) {
      group <- by_source[[source]]
      ancestors <- vapply(group, function(lm) pop_names[lm$dest + 1], "")
      proportions <- .sequential_to_proportions(vapply(group, function(lm) lm$proportion, 0))
      if (.isclose(sum(proportions), 1)) {
        demes[[source]]$ancestors <- as.list(ancestors)
        demes[[source]]$proportions <- as.list(proportions)
        demes[[source]]$start_time <- time
        for (a in ancestors) {
          if (!isTRUE(resolved$get(a)$initially_active)) {
            n_ep <- length(demes[[a]]$epochs)
            demes[[a]]$epochs[[n_ep]]$end_time <- time
          }
        }
      } else {
        pulse_keys <- c(pulse_keys, vapply(group, function(lm) {
          paste(lm$source, lm$dest, format(lm$proportion, digits = 17))
        }, ""))
      }
    }
    for (lm in lineage_movements) {
      if (paste(lm$source, lm$dest, format(lm$proportion, digits = 17)) %in% pulse_keys) {
        pulses[[length(pulses) + 1]] <- list(sources = list(pop_names[lm$dest + 1]),
                                             dest = pop_names[lm$source + 1],
                                             proportions = list(lm$proportion), time = time)
      }
    }
  }
  for (name in names(demes)) {
    deme <- demes[[name]]
    start_time <- deme$start_time %||% Inf
    epoch <- deme$epochs[[1]]
    growth_rate <- epoch$growth_rate %||% 0
    epoch$growth_rate <- NULL
    if (growth_rate != 0) {
      if (is.infinite(start_time)) {
        .stopf("%s: growth rate for infinite-length epoch is invalid", name)
      }
      epoch$start_size <- epoch$end_size * exp(-(start_time - epoch$end_time) * growth_rate)
    } else {
      epoch$start_size <- epoch$end_size
    }
    deme$epochs[[1]] <- epoch
    demes[[name]] <- deme
  }
  mm_list <- lapply(rev(dbg$epochs), function(e) t(e$migration_matrix))
  end_times <- vapply(rev(dbg$epochs), function(e) e$start_time, numeric(1))
  migrations <- list()
  current <- list()
  start_time <- Inf
  n <- length(pop_names)
  for (q in seq_along(mm_list)) {
    mm <- mm_list[[q]]
    end_time <- end_times[q]
    for (j in seq_len(n)) {
      for (k in seq_len(n)) {
        if (j == k) next
        rate <- mm[j, k]
        key <- paste(j, k)
        cur <- current[[key]]
        if (is.null(cur)) {
          if (rate != 0) {
            migrations[[length(migrations) + 1]] <- list(
              source = pop_names[j], dest = pop_names[k], start_time = start_time,
              end_time = end_time, rate = rate)
            current[[key]] <- length(migrations)
          }
        } else if (rate == 0) {
          current[[key]] <- NULL
        } else if (migrations[[cur]]$rate == rate) {
          migrations[[cur]]$end_time <- end_time
        } else {
          migrations[[length(migrations) + 1]] <- list(
            source = pop_names[j], dest = pop_names[k], start_time = start_time,
            end_time = end_time, rate = rate)
          current[[key]] <- length(migrations)
        }
      }
    }
    start_time <- end_time
  }
  starts <- vapply(demes, function(d) d$start_time %||% Inf, numeric(1))
  demes <- unname(demes[order(-starts, method = "radix")])
  data <- list(time_units = "generations", demes = demes, migrations = migrations,
               pulses = rev(pulses))
  if (length(migrations) == 0) data$migrations <- NULL
  if (length(pulses) == 0) data$pulses <- NULL
  resolve_demes(data)
}

#' Read and write demes models
#'
#' `read_demes()` reads a demes YAML file and returns the resolved graph;
#' `write_demes()` writes a graph (e.g. from `demography$to_demes()`) as YAML.
#' Both require the 'yaml' package.
#'
#' @param file A file path.
#' @param graph A demes graph.
#' @return `read_demes()` returns a `demes_graph`; `write_demes()` returns
#'   `file` invisibly.
#' @export
read_demes <- function(file) {
  resolve_demes(.read_yaml(file))
}

#' @rdname read_demes
#' @export
write_demes <- function(graph, file) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The 'yaml' package is required to write demes models", call. = FALSE)
  }
  out <- unclass(graph)
  out$demes <- lapply(out$demes, function(d) {
    d$epochs <- lapply(d$epochs, function(e) {
      e$start_time <- NULL
      e
    })
    if (length(d$ancestors) == 0) {
      d$ancestors <- list()
      d$proportions <- list()
    } else {
      d$ancestors <- as.list(d$ancestors)
      d$proportions <- as.list(d$proportions)
    }
    d
  })
  out$pulses <- lapply(out$pulses, function(p) {
    p$sources <- as.list(p$sources)
    p$proportions <- as.list(p$proportions)
    p
  })
  yaml::write_yaml(out, file)
  invisible(file)
}
