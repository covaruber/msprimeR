# Internal utilities shared across the package.

`%||%` <- function(x, y) if (is.null(x)) y else x

# Python's builtin sum(): sequential double precision accumulation (R's sum()
# accumulates in long double, which can differ in the last bits).
.psum <- function(x) {
  s <- 0
  for (v in x) s <- s + v
  s
}

# numpy's pairwise summation of a contiguous float64 vector.
.np_sum <- function(x) {
  n <- length(x)
  if (n < 8) return(.psum(x))
  if (n > 128) {
    n2 <- (n %/% 2) - ((n %/% 2) %% 8)
    return(.np_sum(x[seq_len(n2)]) + .np_sum(x[(n2 + 1):n]))
  }
  r <- x[1:8]
  i <- 9
  limit <- n - (n %% 8)
  while (i <= limit) {
    r <- r + x[i:(i + 7)]
    i <- i + 8
  }
  res <- ((r[1] + r[2]) + (r[3] + r[4])) + ((r[5] + r[6]) + (r[7] + r[8]))
  while (i <= n) {
    res <- res + x[i]
    i <- i + 1
  }
  res
}

.np_row_sums <- function(m) vapply(seq_len(nrow(m)), function(j) .np_sum(m[j, ]), numeric(1))

.stopf <- function(...) stop(sprintf(...), call. = FALSE)

.warnf <- function(...) warning(sprintf(...), call. = FALSE)

# Parse a boolean flag, which can be TRUE, FALSE or NULL (the default).
.parse_flag <- function(value, default) {
  if (is.null(value)) {
    return(default)
  }
  if (!is.logical(value) || length(value) != 1 || is.na(value)) {
    stop("Boolean flag must be TRUE, FALSE, or NULL (the default value)", call. = FALSE)
  }
  value
}

# TRUE if value can be converted losslessly to an integer.
.isinteger <- function(value) {
  is.numeric(value) && length(value) == 1 && !is.na(value) && is.finite(value) &&
    value == round(value)
}

.is_number <- function(value) {
  is.numeric(value) && length(value) == 1 && !is.na(value)
}

.check_number <- function(value, name) {
  if (!.is_number(value)) .stopf("'%s' must be a single number", name)
  as.double(value)
}

# Random seeds are drawn from R's RNG so that set.seed() makes unseeded
# simulations reproducible. Valid seeds are in [1, 2^32 - 1].
.get_random_seed <- function() {
  floor(stats::runif(1, 1, 2^32))
}

.parse_random_seed <- function(seed) {
  if (is.null(seed)) {
    seed <- .get_random_seed()
  }
  if (!.isinteger(seed)) stop("random_seed must be an integer", call. = FALSE)
  seed <- as.double(seed)
  if (seed < 1 || seed >= 2^32) {
    stop("seeds must be greater than 0 and less than 2^32", call. = FALSE)
  }
  seed
}

# Pack a list/vector of strings into raw bytes + offsets (tskit ragged layout).
.pack_strings <- function(x) {
  x <- as.character(x)
  raws <- lapply(x, function(s) charToRaw(enc2utf8(s)))
  .pack_raw(raws)
}

.pack_raw <- function(raws) {
  lens <- vapply(raws, length, numeric(1))
  list(data = if (length(raws)) unlist(raws, use.names = FALSE) else raw(0),
       offset = c(0, cumsum(lens)))
}

.unpack_strings <- function(data, offset) {
  .Call(msprimer_unpack_strings, as.raw(data), as.double(offset))
}

.unpack_raw <- function(data, offset) {
  n <- length(offset) - 1
  if (n <= 0) return(list())
  lapply(seq_len(n), function(j) {
    if (offset[j + 1] > offset[j]) data[(offset[j] + 1):offset[j + 1]] else raw(0)
  })
}

# ---------------------------------------------------------------------------
# JSON encoding compatible with Python's json module.
# ---------------------------------------------------------------------------

.json_escape_string <- function(s) {
  s <- enc2utf8(s)
  cps <- utf8ToInt(s)
  if (length(cps) == 1 && is.na(cps)) stop("Invalid UTF-8 string", call. = FALSE)
  out <- character(length(cps))
  for (j in seq_along(cps)) {
    cp <- cps[j]
    out[j] <- if (cp == 0x22) {
      "\\\""
    } else if (cp == 0x5c) {
      "\\\\"
    } else if (cp == 0x0a) {
      "\\n"
    } else if (cp == 0x0d) {
      "\\r"
    } else if (cp == 0x09) {
      "\\t"
    } else if (cp == 0x08) {
      "\\b"
    } else if (cp == 0x0c) {
      "\\f"
    } else if (cp < 0x20 || (cp > 0x7e && cp <= 0xffff)) {
      sprintf("\\u%04x", cp)
    } else if (cp > 0xffff) {
      v <- cp - 0x10000
      sprintf("\\u%04x\\u%04x", 0xd800 + bitwShiftR(v, 10), 0xdc00 + bitwAnd(v, 0x3ff))
    } else {
      intToUtf8(cp)
    }
  }
  paste0("\"", paste(out, collapse = ""), "\"")
}

# Shortest decimal representation that round-trips, as Python's repr(float).
.json_number <- function(x, is_int) {
  if (is.na(x)) return(if (is.nan(x)) "NaN" else "null")
  if (is.infinite(x)) return(if (x > 0) "Infinity" else "-Infinity")
  if (is_int) return(format(x, scientific = FALSE, trim = TRUE))
  # R has no separate integer literal in common use, so whole numbers are
  # written as JSON integers (as Python does for ints).
  if (x == round(x) && abs(x) < 1e15) return(format(x, scientific = FALSE, trim = TRUE))
  s <- NULL
  for (d in 15:17) {
    s <- sprintf(paste0("%.", d, "g"), x)
    if (as.double(s) == x) break
  }
  # Python uses a two digit exponent minimum ('1e-05') like C, and adds '.0'
  # to integral values written without an exponent.
  if (!grepl("[.eE]", s)) s <- paste0(s, ".0")
  s
}

.json_encode <- function(x, sort_keys = TRUE, item_sep = ",", key_sep = ":") {
  enc <- function(v) .json_encode(v, sort_keys, item_sep, key_sep)
  if (is.null(x)) return("null")
  if (inherits(x, "R6")) {
    if (is.function(x$asdict)) return(enc(x$asdict()))
    stop("Cannot JSON encode R6 object without an asdict() method", call. = FALSE)
  }
  if (is.environment(x) || is.function(x)) {
    stop("Cannot JSON encode environments or functions", call. = FALSE)
  }
  if (is.list(x)) {
    nms <- names(x)
    if (!is.null(nms) && (length(x) > 0 || inherits(x, "json_object"))) {
      idx <- seq_along(x)
      if (sort_keys) idx <- order(nms, method = "radix")
      items <- vapply(idx, function(j) {
        paste0(.json_escape_string(nms[j]), key_sep, enc(x[[j]]))
      }, character(1))
      return(paste0("{", paste(items, collapse = item_sep), "}"))
    }
    if (inherits(x, "json_object")) return("{}")
    items <- vapply(x, enc, character(1))
    return(paste0("[", paste(items, collapse = item_sep), "]"))
  }
  if (is.matrix(x)) {
    return(enc(lapply(seq_len(nrow(x)), function(j) x[j, ])))
  }
  scalar <- length(x) == 1 && !inherits(x, "AsIs")
  vals <- if (is.character(x)) {
    vapply(x, function(s) if (is.na(s)) "null" else .json_escape_string(s), character(1))
  } else if (is.logical(x)) {
    ifelse(is.na(x), "null", ifelse(x, "true", "false"))
  } else if (is.integer(x)) {
    vapply(x, function(v) .json_number(v, TRUE), character(1))
  } else if (is.numeric(x)) {
    vapply(x, function(v) .json_number(v, FALSE), character(1))
  } else if (is.raw(x)) {
    vapply(as.integer(x), function(v) .json_number(v, TRUE), character(1))
  } else if (is.factor(x)) {
    vapply(as.character(x), .json_escape_string, character(1))
  } else {
    stop(sprintf("Cannot JSON encode object of class '%s'", class(x)[1]), call. = FALSE)
  }
  if (scalar) return(unname(vals))
  paste0("[", paste(vals, collapse = item_sep), "]")
}

#' Canonical JSON encoding
#'
#' Encodes an R object as JSON in the canonical form used by tskit for metadata
#' and metadata schemas (sorted keys, no whitespace, ASCII output), identical to
#' Python's `json.dumps(obj, sort_keys=True, separators=(",", ":"))`.
#' Length-one atomic vectors are encoded as scalars; wrap them in [I()] to
#' force an array.
#'
#' @param x An R object (named lists become objects, unnamed lists arrays).
#' @return A character scalar.
#' @examples
#' canonical_json(list(name = "pop_0", description = ""))
#' @export
canonical_json <- function(x) {
  .json_encode(x, sort_keys = TRUE)
}

.json_decode <- function(s) {
  if (is.raw(s)) s <- rawToChar(s)
  if (!nzchar(s)) return(NULL)
  jsonlite::fromJSON(s, simplifyVector = TRUE, simplifyDataFrame = FALSE,
                     simplifyMatrix = FALSE)
}

.json_object <- function(x = list()) {
  if (length(x) == 0) x <- structure(list(), names = character(0))
  structure(x, class = c("json_object", class(x)))
}
