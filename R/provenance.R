# Provenance records conforming to the tskit provenance schema.

.package_version <- function() {
  as.character(utils::packageVersion("msprimeR"))
}

#' Versions of the bundled C libraries
#'
#' @return A named list with the GSL, tskit and kastore versions.
#' @export
library_versions <- function() {
  v <- .Call(msprimer_versions)
  lapply(v, function(x) paste(x, collapse = "."))
}

.environment_cache <- new.env(parent = emptyenv())

.get_environment <- function() {
  if (is.null(.environment_cache$env)) {
    info <- Sys.info()
    libs <- library_versions()
    .environment_cache$env <- list(
      os = list(system = unname(info[["sysname"]]), node = unname(info[["nodename"]]),
                release = unname(info[["release"]]), version = unname(info[["version"]]),
                machine = unname(info[["machine"]])),
      R = list(implementation = "R", version = paste(R.version$major, R.version$minor,
                                                    sep = ".")),
      libraries = list(gsl = list(version = libs$gsl),
                       tskit = list(version = libs$tskit),
                       kastore = list(version = libs$kastore)))
  }
  .environment_cache$env
}

.get_provenance_dict <- function(parameters) {
  list(
    schema_version = "1.0.0",
    software = list(name = "msprimeR", version = .package_version()),
    parameters = parameters,
    environment = .get_environment())
}

# Convert parameter values into JSON-encodable form, mirroring msprime's
# ProvenanceEncoderDecoder: objects are encoded via their constructor arguments
# along with a "__class__" field; tree sequences are not stored.
.provenance_value <- function(x) {
  if (is.null(x)) return(NULL)
  if (.is_tree_sequence(x) || inherits(x, "TableCollection")) {
    return(list(`__constant__` = "__current_ts__"))
  }
  if (inherits(x, "R6")) {
    if (is.function(x$asdict)) {
      d <- lapply(x$asdict(), .provenance_value)
      d[["__class__"]] <- paste0("msprimeR.", class(x)[1])
      return(d)
    }
    return(list(`__error__` = paste("Cannot encode", class(x)[1])))
  }
  if (is.function(x)) {
    return(list(`__error__` = "Functions are not serialised in provenance"))
  }
  if (inherits(x, "msprimeR_object")) {
    d <- lapply(unclass(x), .provenance_value)
    d[["__class__"]] <- paste0("msprimeR.", class(x)[1])
    return(d)
  }
  if (is.list(x)) {
    out <- lapply(x, .provenance_value)
    if (!is.null(names(x))) out <- .json_object(out)
    return(out)
  }
  if (is.matrix(x)) return(lapply(seq_len(nrow(x)), function(j) x[j, ]))
  if (length(x) != 1 && is.atomic(x)) return(I(x))
  x
}

.json_encode_provenance <- function(provenance_dict) {
  prov <- .provenance_value(provenance_dict)
  .json_encode(prov, sort_keys = FALSE, item_sep = ", ", key_sep = ": ")
}
