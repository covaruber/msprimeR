# Tree sequence results and conversion to/from the tskit data model.
#
# Simulation results are tree sequences stored as a table collection in the
# same "dict" layout used by tskit's TableCollection.asdict() (and tskitR).

.empty_tables <- function(sequence_length) {
  empty_meta <- list(metadata = raw(0), metadata_offset = 0, metadata_schema = "")
  list(
    encoding_version = c(1L, 6L),
    sequence_length = as.double(sequence_length),
    metadata_schema = "",
    metadata = raw(0),
    time_units = "unknown",
    individuals = c(list(flags = numeric(0), location = numeric(0), location_offset = 0,
                         parents = integer(0), parents_offset = 0), empty_meta),
    nodes = c(list(flags = numeric(0), time = numeric(0), population = integer(0),
                   individual = integer(0)), empty_meta),
    edges = c(list(left = numeric(0), right = numeric(0), parent = integer(0),
                   child = integer(0)), empty_meta),
    migrations = c(list(left = numeric(0), right = numeric(0), node = integer(0),
                        source = integer(0), dest = integer(0), time = numeric(0)),
                   empty_meta),
    sites = c(list(position = numeric(0), ancestral_state = raw(0),
                   ancestral_state_offset = 0), empty_meta),
    mutations = c(list(site = integer(0), node = integer(0), time = numeric(0),
                       derived_state = raw(0), derived_state_offset = 0,
                       parent = integer(0)), empty_meta),
    populations = empty_meta,
    provenances = list(timestamp = raw(0), timestamp_offset = 0, record = raw(0),
                       record_offset = 0)
  )
}

.num_rows <- function(table, col) length(table[[col]])

# Append rows to a dict table: `cols` are the fixed columns, `ragged` a named
# list of lists of vectors (one per row) for the ragged columns.
.append_rows <- function(table, cols, ragged = list()) {
  n <- NULL
  for (nm in names(cols)) {
    table[[nm]] <- c(table[[nm]], cols[[nm]])
    n <- length(cols[[nm]])
  }
  for (nm in names(ragged)) {
    rows <- ragged[[nm]]
    off_name <- paste0(nm, "_offset")
    old <- table[[off_name]]
    lens <- vapply(rows, length, numeric(1))
    if (length(rows) > 0) table[[nm]] <- c(table[[nm]], unlist(rows, use.names = FALSE))
    table[[off_name]] <- c(old, old[length(old)] + cumsum(lens))
    n <- length(rows)
  }
  # Ragged columns that were not supplied get empty entries.
  ragged_cols <- sub("_offset$", "", grep("_offset$", names(table), value = TRUE))
  for (nm in setdiff(ragged_cols, names(ragged))) {
    off_name <- paste0(nm, "_offset")
    old <- table[[off_name]]
    table[[off_name]] <- c(old, rep(old[length(old)], n))
  }
  table
}

.add_provenance <- function(tables, record) {
  # Same format as Python's datetime.now().isoformat() used by tskit.
  timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6")
  tables$provenances <- .append_rows(
    tables$provenances, list(),
    list(timestamp = list(charToRaw(timestamp)), record = list(charToRaw(enc2utf8(record)))))
  tables
}

.new_ts <- function(tables) {
  res <- .Call(msprimer_ts_summary, tables)
  structure(list(tables = res$tables, summary = res$summary), class = "msprime_ts")
}

#' Tree sequence objects
#'
#' Simulation functions return tree sequences. By default (option
#' `msprimeR.output = "auto"`) these are `tskitR::TreeSequence` objects when the
#' 'tskitR' package is installed, giving access to the full 'tskit' API
#' (trees, statistics, genotypes, ...). Otherwise, or when
#' `options(msprimeR.output = "msprime_ts")` is set, a lightweight
#' `msprime_ts` object is returned. This holds the table collection in the
#' standard tskit layout (`ts$tables`), basic properties (e.g. `ts$num_trees`,
#' `ts$num_nodes`) and can be written to the native `.trees` format with
#' [ts_dump()] or converted with [as_tskit()].
#'
#' All functions in msprimeR accepting a tree sequence work with either kind of
#' object, with a `tskitR::TableCollection`, or with the path of a `.trees`
#' file.
#'
#' @name msprime_ts
NULL

.ts_fields <- c("num_trees", "num_samples", "num_nodes", "num_edges", "num_sites",
                "num_mutations", "num_individuals", "num_populations",
                "num_migrations", "num_provenances", "sequence_length",
                "discrete_genome", "discrete_time", "min_time", "max_time", "time_units")

#' @export
`$.msprime_ts` <- function(x, name) {
  if (name %in% c("tables", "summary")) return(.subset2(x, name))
  if (name %in% .ts_fields) return(.subset2(x, "summary")[[name]])
  NULL
}

#' @export
print.msprime_ts <- function(x, ...) {
  s <- .subset2(x, "summary")
  cat("<msprime_ts> tree sequence\n")
  rows <- c("num_trees", "sequence_length", "time_units", "num_samples", "num_nodes",
            "num_edges", "num_sites", "num_mutations", "num_individuals",
            "num_populations", "num_migrations", "num_provenances")
  w <- max(nchar(rows))
  for (r in rows) cat(sprintf("  %-*s  %s\n", w, r, format(s[[r]])))
  invisible(x)
}

.output_mode <- function() {
  mode <- getOption("msprimeR.output", "auto")
  if (!mode %in% c("auto", "tskitR", "msprime_ts")) {
    stop("option msprimeR.output must be one of 'auto', 'tskitR' or 'msprime_ts'",
         call. = FALSE)
  }
  if (mode == "auto") {
    mode <- if (requireNamespace("tskitR", quietly = TRUE)) "tskitR" else "msprime_ts"
  }
  mode
}

# Build the user-facing result from a table dict.
.make_result <- function(tables) {
  ts <- .new_ts(tables)
  if (.output_mode() == "tskitR") as_tskit(ts) else ts
}

#' Convert to a tskitR tree sequence
#'
#' Converts an `msprime_ts` (or any object accepted by [ts_tables()]) into a
#' `tskitR::TreeSequence`. Requires the 'tskitR' package.
#'
#' @param ts A tree sequence.
#' @return A `tskitR::TreeSequence`.
#' @export
as_tskit <- function(ts) {
  if (!requireNamespace("tskitR", quietly = TRUE)) {
    stop("The 'tskitR' package is required for as_tskit()", call. = FALSE)
  }
  if (inherits(ts, "TreeSequence")) return(ts)
  tables <- ts_tables(ts)
  tskitR::table_collection_fromdict(tables)$tree_sequence()
}

#' Convert to an msprime_ts object
#'
#' @param ts A tree sequence (`tskitR::TreeSequence`, `msprime_ts`, table
#'   collection or `.trees` file path).
#' @return An `msprime_ts` object.
#' @export
as_msprime_ts <- function(ts) {
  if (inherits(ts, "msprime_ts")) return(ts)
  .new_ts(ts_tables(ts))
}

#' Table collection of a tree sequence
#'
#' Returns the tables of a tree sequence as a named list in the layout of
#' tskit's `TableCollection.asdict()`.
#'
#' @param ts An `msprime_ts`, `tskitR::TreeSequence`, `tskitR::TableCollection`,
#'   a table dict, or the path of a `.trees` file.
#' @return A named list.
#' @export
ts_tables <- function(ts) {
  if (inherits(ts, "msprime_ts")) return(.subset2(ts, "tables"))
  if (is.character(ts) && length(ts) == 1) return(.Call(msprimer_tables_load, ts))
  if (inherits(ts, "R6")) {
    if (is.function(ts$dump_tables)) return(.fix_tskitr_dict(ts$dump_tables()$asdict()))
    if (is.function(ts$asdict)) return(.fix_tskitr_dict(ts$asdict()))
  }
  if (is.list(ts) && !is.null(ts$sequence_length) && !is.null(ts$nodes)) {
    return(.Call(msprimer_tables_roundtrip, ts))
  }
  stop("Expected a tree sequence: an msprime_ts, a tskitR TreeSequence or ",
       "TableCollection, a table dict, or the path of a .trees file", call. = FALSE)
}

# tskitR asdict() output is already in the right layout; normalise it in C.
.fix_tskitr_dict <- function(d) {
  .Call(msprimer_tables_roundtrip, d)
}

.is_tree_sequence <- function(x) {
  inherits(x, "msprime_ts") || inherits(x, "TreeSequence")
}

#' Write a tree sequence to a .trees file
#'
#' Writes the tree sequence in the native 'kastore' based `.trees` format, which
#' can be read by 'tskitR' (`tskitR::tskit_load()`), Python's
#' `tskit.load()` and any other tskit-compatible software.
#'
#' @param ts A tree sequence.
#' @param file The file path.
#' @return `file`, invisibly.
#' @export
ts_dump <- function(ts, file) {
  .Call(msprimer_tables_dump, ts_tables(ts), path.expand(file))
  invisible(file)
}

#' Read a tree sequence from a .trees file
#'
#' @param file The file path.
#' @return A tree sequence (see [msprime_ts] for the class returned).
#' @export
ts_load <- function(file) {
  .make_result(.Call(msprimer_tables_load, path.expand(file)))
}

#' Extract a table as a data frame
#'
#' @param ts A tree sequence (or table dict).
#' @param name One of "individuals", "nodes", "edges", "migrations", "sites",
#'   "mutations", "populations" or "provenances".
#' @return A data frame; text columns are decoded to character, metadata is
#'   returned as a list column of raw vectors.
#' @export
ts_table <- function(ts, name) {
  tables <- ts_tables(ts)
  t <- tables[[name]]
  if (is.null(t)) .stopf("Unknown table '%s'", name)
  out <- list()
  offsets <- grep("_offset$", names(t), value = TRUE)
  ragged <- sub("_offset$", "", offsets)
  n <- length(t[[offsets[1]]]) - 1
  for (nm in names(t)) {
    if (nm %in% offsets || nm == "metadata_schema") next
    if (nm %in% ragged) {
      v <- t[[nm]]
      off <- t[[paste0(nm, "_offset")]]
      if (nm %in% c("ancestral_state", "derived_state", "timestamp", "record")) {
        out[[nm]] <- .unpack_strings(v, off)
      } else if (nm == "metadata") {
        out[[nm]] <- I(.unpack_raw(v, off))
      } else {
        out[[nm]] <- I(lapply(seq_len(n), function(j) {
          if (off[j + 1] > off[j]) v[(off[j] + 1):off[j + 1]] else v[0]
        }))
      }
    } else {
      out[[nm]] <- t[[nm]]
    }
  }
  df <- as.data.frame(out, stringsAsFactors = FALSE)
  if (nrow(df) != n) df <- df[seq_len(n), , drop = FALSE]
  df
}
