# Rate maps: piecewise-constant rates along the genome (port of
# msprime/intervals.py).

#' RateMap: a piecewise-constant rate along the genome
#'
#' A `RateMap` maps a non-negative rate to each of a set of contiguous
#' intervals `[position[j], position[j + 1])` covering the genome. Intervals
#' with unknown rate (missing data) are encoded by `NaN`.
#'
#' Create rate maps with `RateMap$new(position, rate)`,
#' `RateMap$uniform(sequence_length, rate)` or
#' `RateMap$read_hapmap(file, ...)`.
#'
#' Interval indices (e.g. returned by `find_index()`) are 1-based, so that
#' `rate_map$rate[rate_map$find_index(x)]` gives the rate at `x`.
#'
#' @field position The `n + 1` breakpoints, starting at 0 and ending at the
#'   sequence length.
#' @field rate The `n` rates (NaN for missing data).
#' @field left,right,mid,span Left/right coordinates, midpoints and spans of the
#'   intervals.
#' @field mass Rate times span for each interval.
#' @field missing,non_missing Logical vectors flagging missing intervals.
#' @field num_intervals,num_missing_intervals,num_non_missing_intervals Counts.
#' @field sequence_length The sequence length covered by the map.
#' @field total_mass The total mass of the map.
#' @field mean_rate The span-weighted mean rate over non-missing intervals.
#' @export
RateMap <- R6::R6Class("RateMap",
  cloneable = FALSE,
  public = list(
    #' @description Create a new rate map.
    #' @param position A vector of `n + 1` positions, starting at 0.
    #' @param rate A vector of `n` non-negative rates (NaN for missing data).
    initialize = function(position, rate) {
      position <- as.double(position)
      rate <- as.double(rate)
      size <- length(position)
      if (size < 2) stop("Must have at least two positions", call. = FALSE)
      if (length(rate) != size - 1) {
        stop("Rate array must have one less entry than the position array", call. = FALSE)
      }
      if (any(is.na(position))) stop("Positions must not be missing", call. = FALSE)
      if (position[1] != 0) stop("First position must be zero", call. = FALSE)
      span <- diff(position)
      if (any(span <= 0)) {
        .stopf("Position values not strictly increasing at indexes %s",
               paste(which(span <= 0) + 1, collapse = ", "))
      }
      if (any(rate < 0, na.rm = TRUE)) {
        .stopf("Rate values negative at indexes %s",
               paste(which(rate < 0), collapse = ", "))
      }
      private$.position <- position
      private$.rate <- rate
      private$.missing <- is.nan(rate) | is.na(rate)
      if (all(private$.missing)) stop("All intervals are missing data", call. = FALSE)
      mass <- rate * span
      mass[private$.missing] <- 0
      private$.cumulative_mass <- c(0, cumsum(mass))
    },
    #' @description Rates at the specified positions (NaN in missing intervals).
    #' @param x Numeric positions.
    get_rate = function(x) {
      loc <- findInterval(x, private$.position)
      if (any(loc < 1) || any(loc > length(private$.rate))) {
        stop("position out of bounds", call. = FALSE)
      }
      private$.rate[loc]
    },
    #' @description Cumulative mass up to (but not including) each position.
    #' @param x Numeric positions.
    get_cumulative_mass = function(x) {
      x <- as.double(x)
      if (any(x < 0) || any(x > self$sequence_length)) {
        .stopf("Cannot have positions < 0 or > %s", format(self$sequence_length))
      }
      stats::approx(private$.position, private$.cumulative_mass, xout = x,
                    rule = 2, ties = "ordered")$y
    },
    #' @description The (1-based) index of the interval containing `x`.
    #' @param x A single position.
    find_index = function(x) {
      if (x < 0 || x >= self$sequence_length) .stopf("Position %s out of bounds", x)
      findInterval(x, private$.position)
    },
    #' @description The coordinates of missing intervals, as a two column matrix.
    missing_intervals = function() {
      cbind(left = self$left[private$.missing], right = self$right[private$.missing])
    },
    #' @description List representation (position, rate).
    asdict = function() list(position = private$.position, rate = private$.rate),
    #' @description A copy of this map.
    copy = function() RateMap$new(private$.position, private$.rate),
    #' @description A subset of this map in `[left, right)`.
    #' @param left Left coordinate (default 0).
    #' @param right Right coordinate (default sequence length).
    #' @param trim If TRUE, shift coordinates so the new map starts at zero; if
    #'   FALSE, keep the coordinate system and mark flanks as missing.
    slice = function(left = NULL, right = NULL, trim = FALSE) {
      L <- self$sequence_length
      left <- left %||% 0
      right <- right %||% L
      if (!(0 <= left && left < right && right <= L)) {
        .stopf("Invalid slice: left=%s, right=%s", left, right)
      }
      pos <- private$.position
      i <- self$find_index(left)
      j <- i - 1 + findInterval(right, pos[i:length(pos)], left.open = FALSE)
      # j is now the number of positions <= right (Python searchsorted 'right').
      if (right > pos[j]) j <- j + 1
      position <- pos[i:j]
      rate <- private$.rate[i:(j - 1)]
      position[1] <- left
      position[length(position)] <- right
      if (trim) return(RateMap$new(position - left, rate))
      if (left != 0) {
        if (is.na(rate[1])) {
          position[1] <- 0
        } else {
          rate <- c(NaN, rate)
          position <- c(0, position)
        }
      }
      if (right != L) {
        if (is.na(rate[length(rate)])) {
          position[length(position)] <- L
        } else {
          rate <- c(rate, NaN)
          position <- c(position, L)
        }
      }
      RateMap$new(position, rate)
    },
    #' @description Print the map.
    #' @param ... Ignored.
    print = function(...) {
      n <- self$num_intervals
      idx <- if (n < 40) seq_len(n) else c(1:10, (n - 9):n)
      df <- data.frame(left = self$left[idx], right = self$right[idx],
                       mid = self$mid[idx], span = self$span[idx],
                       rate = self$rate[idx])
      cat(sprintf("<RateMap> %d intervals, sequence_length = %s\n", n,
                  format(self$sequence_length)))
      if (n >= 40) {
        print(utils::head(df, 10), row.names = FALSE)
        cat("  ...\n")
        print(utils::tail(df, 10), row.names = FALSE)
      } else {
        print(df, row.names = FALSE)
      }
      invisible(self)
    }
  ),
  private = list(.position = NULL, .rate = NULL, .missing = NULL,
                 .cumulative_mass = NULL),
  active = list(
    position = function() private$.position,
    rate = function() private$.rate,
    left = function() private$.position[-length(private$.position)],
    right = function() private$.position[-1],
    mid = function() self$left + self$span / 2,
    span = function() diff(private$.position),
    mass = function() private$.rate * self$span,
    missing = function() private$.missing,
    non_missing = function() !private$.missing,
    num_intervals = function() length(private$.rate),
    num_missing_intervals = function() sum(private$.missing),
    num_non_missing_intervals = function() sum(!private$.missing),
    sequence_length = function() private$.position[length(private$.position)],
    total_mass = function() private$.cumulative_mass[length(private$.cumulative_mass)],
    mean_rate = function() self$total_mass / sum(self$span[!private$.missing])
  )
)

RateMap$uniform <- function(sequence_length, rate) {
  RateMap$new(position = c(0, sequence_length), rate = rate)
}

RateMap$read_hapmap <- function(file, sequence_length = NULL, has_header = TRUE,
                                position_col = NULL, rate_col = NULL, map_col = NULL) {
  read_hapmap(file, sequence_length = sequence_length, has_header = has_header,
              position_col = position_col, rate_col = rate_col, map_col = map_col)
}

#' Read a HapMap format recombination map
#'
#' Parses a whitespace-delimited file in HapMap format and returns a
#' [RateMap]. By default the second column (zero-based index 1) is the
#' physical position in base pairs and the fourth column (index 3) the genetic
#' map position in centiMorgans.
#'
#' If the first physical position is greater than zero and the first genetic
#' position is zero, the region before it is marked as missing (NaN). If
#' `sequence_length` is larger than the last physical position, the trailing
#' region is also marked as missing.
#'
#' @param file The file to read (may be compressed).
#' @param sequence_length Optional total length of the map.
#' @param has_header If TRUE (default) the first line is skipped.
#' @param position_col Zero-based column index of the physical position.
#' @param rate_col Zero-based column index of the rate in cM/Mb (if used, the
#'   last value in this column must be zero).
#' @param map_col Zero-based column index of the genetic position in cM.
#' @return A [RateMap].
#' @export
read_hapmap <- function(file, sequence_length = NULL, has_header = TRUE,
                        position_col = NULL, rate_col = NULL, map_col = NULL) {
  if (is.null(rate_col) && is.null(map_col)) {
    map_col <- 3
  } else if (!is.null(rate_col) && !is.null(map_col)) {
    stop("Cannot specify both rate_col and map_col", call. = FALSE)
  }
  position_col <- position_col %||% 1
  value_col <- if (!is.null(map_col)) map_col else rate_col
  if (position_col == value_col) {
    stop("Cannot specify the same columns for position_col and rate_col or map_col",
         call. = FALSE)
  }
  data <- utils::read.table(file, header = FALSE, skip = if (has_header) 1 else 0,
                            colClasses = "character", fill = TRUE)
  if (nrow(data) == 0) stop("Empty hapmap file", call. = FALSE)
  pos <- as.double(data[[position_col + 1]])
  if (any(pos != round(pos))) stop("Physical positions must be integers", call. = FALSE)
  if (!is.null(map_col)) {
    genetic <- as.double(data[[map_col + 1]]) / 100
  } else {
    rate <- as.double(data[[rate_col + 1]])
    if (rate[length(rate)] != 0) {
      stop("The last entry in the 'rate' column must be zero", call. = FALSE)
    }
    pos_mb <- pos / 1e6
    genetic <- c(0, cumsum(rate[-length(rate)] * diff(pos_mb))) / 100
  }
  start <- pos[1]
  end <- pos[length(pos)]
  if (genetic[1] > 0 && start == 0) {
    stop("The map distance at the start of the chromosome must be zero", call. = FALSE)
  }
  if (start > 0) {
    pos <- c(0, pos)
    if (genetic[1] > 0) start <- 0
    genetic <- c(0, genetic)
  }
  if (!is.null(sequence_length)) {
    if (sequence_length < end) {
      .stopf("The sequence_length cannot be less than the last physical position (%s)",
             format(pos[length(pos)]))
    }
    if (sequence_length > end) {
      pos <- c(pos, sequence_length)
      genetic <- c(genetic, genetic[length(genetic)])
    }
  }
  rate <- diff(genetic) / diff(pos)
  if (start != 0) rate[1] <- NaN
  if (end != pos[length(pos)]) rate[length(rate)] <- NaN
  RateMap$new(position = pos, rate = rate)
}

# Parse a rate parameter (number, NULL or RateMap) into a RateMap.
.parse_rate_map <- function(rate_param, sequence_length, name) {
  msg_head <- sprintf("Error in parsing rate map for %s: ", name)
  if (inherits(rate_param, "RateMap")) {
    if (rate_param$sequence_length != sequence_length) {
      stop(msg_head, "sequence_length must match", call. = FALSE)
    }
    return(rate_param)
  }
  rate <- if (is.null(rate_param)) 0 else .check_number(rate_param, name)
  RateMap$uniform(sequence_length, rate)
}
