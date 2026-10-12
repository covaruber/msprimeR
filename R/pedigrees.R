# Pedigree utilities (port of msprime/pedigrees.py).

# Appends a block of individuals (two nodes each) to a table dict in one go.
.append_pedigree_rows <- function(tables, parents, time, is_sample, population, metadata) {
  N <- nrow(parents)
  if (N == 0) return(list(tables = tables, ids = integer(0)))
  first <- length(tables$individuals$flags)
  ind_ids <- first + seq_len(N) - 1L
  time <- rep_len(as.double(time), N)
  flags <- ifelse(rep_len(is_sample, N), 1, 0)
  population <- rep_len(as.integer(population), N)
  ind <- tables$individuals
  ind$flags <- c(ind$flags, rep(0, N))
  ind$parents <- c(ind$parents, as.integer(t(parents)))
  ind$parents_offset <- c(ind$parents_offset, ind$parents_offset[length(ind$parents_offset)] +
                            2 * seq_len(N))
  ind$location_offset <- c(ind$location_offset,
                           rep(ind$location_offset[length(ind$location_offset)], N))
  md_lens <- if (is.null(metadata)) rep(0, N) else vapply(metadata, length, numeric(1))
  if (!is.null(metadata)) ind$metadata <- c(ind$metadata, unlist(metadata, use.names = FALSE))
  ind$metadata_offset <- c(ind$metadata_offset,
                           ind$metadata_offset[length(ind$metadata_offset)] + cumsum(md_lens))
  tables$individuals <- ind
  tables$nodes <- .append_rows(tables$nodes, list(
    flags = rep(flags, each = 2), time = rep(time, each = 2),
    population = rep(population, each = 2), individual = rep(as.integer(ind_ids), each = 2)))
  list(tables = tables, ids = ind_ids)
}

#' Pedigree builder
#'
#' Builds pedigrees in the format required by the [FixedPedigree()] ancestry
#' model. Individual IDs are zero-based.
#'
#' @field demography The [Demography] defining the populations.
#' @examples
#' pb <- PedigreeBuilder$new()
#' mom <- pb$add_individual(time = 1)
#' dad <- pb$add_individual(time = 1)
#' pb$add_individual(time = 0, parents = c(mom, dad), is_sample = TRUE)
#' tables <- pb$finalise(sequence_length = 100)
#' @export
PedigreeBuilder <- R6::R6Class("PedigreeBuilder",
  cloneable = FALSE,
  public = list(
    demography = NULL,
    #' @description Create a builder.
    #' @param demography A [Demography] (default: one population).
    #' @param individuals_metadata_schema A metadata schema (named list or JSON
    #'   string) for the individuals table; only JSON codecs are supported.
    initialize = function(demography = NULL, individuals_metadata_schema = NULL) {
      demography <- demography %||% Demography$isolated_model(1)
      self$demography <- demography
      tables <- .empty_tables(0)
      tables$time_units <- "generations"
      tables <- demography$insert_populations(tables)
      if (!is.null(individuals_metadata_schema)) {
        tables$individuals$metadata_schema <- if (is.character(individuals_metadata_schema)) {
          individuals_metadata_schema
        } else {
          canonical_json(individuals_metadata_schema)
        }
      }
      private$.tables <- tables
      private$pending <- vector("list", 64)
    },
    #' @description Add an individual, returning its (zero-based) ID.
    #' @param time The time of the individual (generations ago).
    #' @param is_sample Whether the individual's nodes are samples (default:
    #'   `time == 0`).
    #' @param parents Two parent IDs (default `c(-1, -1)`, a founder). Parents
    #'   may refer to individuals not yet added.
    #' @param population The population (name or ID).
    #' @param metadata Metadata: a named list (with a JSON schema) or raw bytes.
    add_individual = function(time, is_sample = NULL, parents = NULL, population = NULL,
                              metadata = NULL) {
      is_sample <- is_sample %||% (time == 0)
      parents <- parents %||% c(-1, -1)
      if (length(parents) != 2) stop("Must have exactly two parents", call. = FALSE)
      md <- private$encode_metadata(metadata)
      population <- private$resolve_population(population)
      if (private$n_pending == length(private$pending)) {
        length(private$pending) <- 2 * length(private$pending)
      }
      private$n_pending <- private$n_pending + 1
      private$pending[[private$n_pending]] <- list(parents = parents, time = time,
                                                   is_sample = is_sample,
                                                   population = population, metadata = md)
      private$num_individuals <- private$num_individuals + 1
      private$num_individuals - 1
    },
    #' @description Add many individuals at once.
    #' @param parents A two column matrix of parent IDs.
    #' @param time The time(s) of the individuals.
    #' @param is_sample Sample status (default: `time == 0`).
    #' @param population The population(s) (zero-based IDs, default 0).
    #' @param metadata NULL or a list of raw vectors (one per individual).
    #' @return The IDs of the new individuals.
    add_individuals = function(parents, time, is_sample = NULL, population = 0,
                               metadata = NULL) {
      private$flush()
      parents <- matrix(as.integer(parents), ncol = 2)
      is_sample <- is_sample %||% (time == 0)
      res <- .append_pedigree_rows(private$.tables, parents, time, is_sample, population,
                                   metadata)
      private$.tables <- res$tables
      private$num_individuals <- private$num_individuals + nrow(parents)
      res$ids
    },
    #' @description The finished pedigree table collection (dict).
    #' @param sequence_length The sequence length (default -1, to be set later).
    finalise = function(sequence_length = NULL) {
      private$flush()
      tables <- private$.tables
      tables$sequence_length <- if (is.null(sequence_length)) -1 else sequence_length
      tables$indexes <- list(edge_insertion_order = integer(0),
                             edge_removal_order = integer(0))
      tables
    }
  ),
  active = list(
    #' @field tables The table collection (dict) built so far.
    tables = function() {
      private$flush()
      private$.tables
    }
  ),
  private = list(
    .tables = NULL,
    pending = NULL,
    n_pending = 0,
    num_individuals = 0,
    flush = function() {
      if (private$n_pending == 0) return(invisible(NULL))
      rows <- private$pending[seq_len(private$n_pending)]
      parents <- matrix(as.integer(unlist(lapply(rows, `[[`, "parents"))), ncol = 2,
                        byrow = TRUE)
      res <- .append_pedigree_rows(
        private$.tables, parents,
        time = vapply(rows, function(r) as.double(r$time), numeric(1)),
        is_sample = vapply(rows, function(r) isTRUE(r$is_sample), logical(1)),
        population = vapply(rows, function(r) as.integer(r$population), integer(1)),
        metadata = lapply(rows, `[[`, "metadata"))
      private$.tables <- res$tables
      private$pending <- vector("list", 64)
      private$n_pending <- 0
      invisible(NULL)
    },
    resolve_population = function(population) {
      if (is.null(population)) {
        if (self$demography$num_populations == 1) return(0L)
        stop("Must specify the population for each individual in a multi-population ",
             "demography", call. = FALSE)
      }
      self$demography$get(population)$id
    },
    encode_metadata = function(metadata) {
      if (is.null(metadata)) return(raw(0))
      if (is.raw(metadata)) return(metadata)
      schema <- private$.tables$individuals$metadata_schema
      if (!nzchar(schema)) {
        stop("Metadata must be raw bytes when no metadata schema is set", call. = FALSE)
      }
      if (!identical(.json_decode(schema)$codec, "json")) {
        stop("Only JSON metadata schemas are supported", call. = FALSE)
      }
      charToRaw(canonical_json(metadata))
    }
  )
)

.parse_pedigree_header <- function(header) {
  tokens <- strsplit(trimws(header), "[[:space:]]+")[[1]]
  if (length(tokens) == 0 || tokens[1] != "#") {
    stop("First line must be the header and start with #", call. = FALSE)
  }
  required <- c("id", "parent0", "parent1", "time")
  tokens <- tokens[-1]
  if (!identical(tokens[seq_len(min(4, length(tokens)))], required)) {
    .stopf("The ['%s'] columns are required", paste(required, collapse = "', '"))
  }
  optional <- tokens[-(1:4)]
  bad <- setdiff(optional, c("is_sample", "population"))
  if (length(bad) > 0) .stopf("Column '%s' not supported", bad[1])
  tokens
}

#' Parse a pedigree file
#'
#' Parses a whitespace-delimited text pedigree with header
#' `# id parent0 parent1 time [is_sample] [population]` into a table
#' collection suitable for the [FixedPedigree()] model. Missing parents are
#' denoted by ".". Individuals receive metadata `{"file_id": <id>}` and are
#' added in file order (parents may appear after their children).
#'
#' @param text_file A file path or connection.
#' @param demography A [Demography] for population references (default: one
#'   population).
#' @param sequence_length If specified, the sequence length of the tables.
#' @return A table collection (named list).
#' @export
parse_pedigree <- function(text_file, demography = NULL, sequence_length = NULL) {
  lines <- readLines(text_file, warn = FALSE)
  if (length(lines) == 0) stop("Pedigree file must contain at least a header", call. = FALSE)
  columns <- .parse_pedigree_header(lines[1])
  rows <- lines[-1]
  line_no <- seq_along(rows) + 1
  keep <- nzchar(trimws(rows))
  rows <- rows[keep]
  line_no <- line_no[keep]
  tokens <- strsplit(trimws(rows), "[[:space:]]+")
  ntok <- lengths(tokens)
  bad <- which(ntok != length(columns))
  if (length(bad) > 0) {
    k <- bad[1]
    .stopf("Incorrect number of columns (found=%d need=%d) at line %d: %s", ntok[k],
           length(columns), line_no[k], rows[k])
  }
  n <- length(rows)
  m <- matrix(unlist(tokens, use.names = FALSE), nrow = n, byrow = TRUE)
  colnames(m) <- columns
  ids <- m[, "id"]
  time <- suppressWarnings(as.double(m[, "time"]))
  if (anyNA(time)) {
    k <- which(is.na(time))[1]
    .stopf("Error parsing value for column 'time' on line %d: could not convert '%s'",
           line_no[k], m[k, "time"])
  }
  is_sample <- NULL
  if ("is_sample" %in% columns) {
    v <- m[, "is_sample"]
    if (any(!v %in% c("0", "1"))) {
      k <- which(!v %in% c("0", "1"))[1]
      .stopf("Error parsing value for column 'is_sample' on line %d: Only 0 or 1 values are supported",
             line_no[k])
    }
    is_sample <- v == "1"
  }
  if (any(ids == ".")) {
    .stopf("'.' cannot be used as an individual ID (line %d)", line_no[which(ids == ".")[1]])
  }
  if (anyDuplicated(ids)) .stopf("Duplicate ID at line %d", line_no[anyDuplicated(ids)])
  builder <- PedigreeBuilder$new(individuals_metadata_schema = list(codec = "json"),
                                 demography = demography)
  parent_ids <- function(col) {
    p <- m[, col]
    idx <- match(p, ids) - 1L
    idx[p == "."] <- -1L
    if (anyNA(idx)) {
      k <- which(is.na(idx))[1]
      .stopf("Parent ID '%s' not defined (line %d)", p[k], line_no[k])
    }
    idx
  }
  parents <- cbind(parent_ids("parent0"), parent_ids("parent1"))
  population <- if ("population" %in% columns) {
    vapply(m[, "population"], function(v) {
      iv <- suppressWarnings(as.integer(v))
      ref <- if (!is.na(iv) && grepl("^[+-]?[0-9]+$", v)) iv else v
      as.integer(builder$demography$get(ref)$id)
    }, integer(1), USE.NAMES = FALSE)
  } else {
    if (builder$demography$num_populations != 1) {
      stop("Must specify the population for each individual in a multi-population ",
           "demography", call. = FALSE)
    }
    0L
  }
  metadata <- lapply(ids, function(id) charToRaw(canonical_json(list(file_id = id))))
  builder$add_individuals(parents = parents, time = time,
                          is_sample = is_sample %||% (time == 0),
                          population = population, metadata = metadata)
  builder$finalise(sequence_length)
}

#' Write a pedigree file
#'
#' Writes the individuals of a pedigree table collection (or tree sequence) in
#' the text format read by [parse_pedigree()].
#'
#' @param ts A tree sequence or table collection.
#' @param out A file path or connection.
#' @return `out`, invisibly.
#' @export
write_pedigree <- function(ts, out) {
  tables <- ts_tables(ts)
  nodes <- tables$nodes
  ind <- tables$individuals
  n_ind <- length(ind$flags)
  by_ind <- split(seq_along(nodes$individual), factor(nodes$individual, levels = seq_len(n_ind) - 1))
  if (any(lengths(by_ind) != 2)) {
    stop("Invalid pedigree format: each individual must be associated with two nodes",
         call. = FALSE)
  }
  first <- vapply(by_ind, `[`, integer(1), 1)
  second <- vapply(by_ind, `[`, integer(1), 2)
  if (any(nodes$time[first] != nodes$time[second])) {
    stop("Pedigree individuals must have the same node time", call. = FALSE)
  }
  s1 <- nodes$flags[first] %% 2
  if (any(s1 != nodes$flags[second] %% 2)) {
    stop("Pedigree individuals must have the same sample status", call. = FALSE)
  }
  if (any(nodes$population[first] != nodes$population[second])) {
    stop("Pedigree individuals must have the same node population", call. = FALSE)
  }
  if (any(diff(ind$parents_offset) != 2)) {
    stop("Each individual must have two parents", call. = FALSE)
  }
  par <- matrix(ind$parents, ncol = 2, byrow = TRUE)
  par <- ifelse(par == -1, ".", as.character(par))
  # is_sample is written as 0/1 (msprime writes True/False, which its own
  # parser rejects).
  time <- vapply(nodes$time[first], .json_number, character(1), is_int = FALSE)
  lines <- c("# id\tparent0\tparent1\ttime\tis_sample\tpopulation",
             sprintf("%d\t%s\t%s\t%s\t%d\t%d", seq_len(n_ind) - 1L, par[, 1], par[, 2],
                     time, as.integer(s1), nodes$population[first]))
  writeLines(lines, out)
  invisible(out)
}

#' Simulate a pedigree
#'
#' Internal utility to simulate a random Wright-Fisher pedigree (as in
#' msprime's `sim_pedigree`). Uses R's random number generator.
#'
#' @param population_size A single size or a function of the generation.
#' @param num_samples Number of samples (default: size at generation 0).
#' @param sequence_length The sequence length of the result.
#' @param random_seed If given, used with [set.seed()] (the global RNG state is
#'   restored afterwards).
#' @param end_time The number of generations.
#' @param direction "forward" or "backward".
#' @return A pedigree table collection.
#' @export
sim_pedigree <- function(population_size = NULL, num_samples = NULL, sequence_length = NULL,
                         random_seed = NULL, end_time = NULL, direction = "forward") {
  N_function <- if (is.function(population_size)) population_size else
    function(t) population_size
  num_samples <- num_samples %||% N_function(0)
  if (!is.null(random_seed)) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old <- if (had_seed) get(".Random.seed", envir = globalenv())
    on.exit({
      if (had_seed) assign(".Random.seed", old, envir = globalenv()) else
        rm(".Random.seed", envir = globalenv())
    })
    set.seed(random_seed)
  }
  choice <- function(population, n) {
    v <- if (length(population) == 1) rep(population, 2 * n) else
      sample(population, 2 * n, replace = TRUE)
    matrix(v, ncol = 2, byrow = TRUE)
  }
  builder <- PedigreeBuilder$new()
  if (direction == "forward") {
    if (num_samples != N_function(0)) {
      stop("if at all specified, num_samples must be equal to population_size at ",
           "generation 0 for forward simulation", call. = FALSE)
    }
    population <- -1L
    for (time in rev(0:end_time)) {
      parents <- choice(population, N_function(time))
      population <- builder$add_individuals(parents = parents, time = time,
                                            is_sample = time == 0)
    }
  } else if (direction == "backward") {
    ancestors <- seq_len(num_samples) - 1L
    parent_offset <- num_samples
    for (time in seq_len(end_time) - 1) {
      population <- seq_len(N_function(time + 1)) - 1L
      parents <- choice(population, length(ancestors))
      unique_parents <- sort(unique(as.vector(parents)))
      parent_ids <- matrix(match(parents, unique_parents) - 1L, ncol = 2)
      builder$add_individuals(parents = parent_ids + parent_offset, time = time,
                              is_sample = time == 0)
      ancestors <- seq_along(unique_parents) - 1L + parent_offset
      parent_offset <- parent_offset + length(ancestors)
    }
    builder$add_individuals(parents = matrix(-1L, length(ancestors), 2), time = end_time,
                            is_sample = end_time == 0)
  } else {
    stop("unknown time direction; choose 'backward' or 'forward'", call. = FALSE)
  }
  builder$finalise(if (is.null(sequence_length)) -1 else sequence_length)
}
