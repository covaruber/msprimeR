# DemographyDebugger (port of msprime.demography.DemographyDebugger).

.POP_STATE <- c("INACTIVE", "ACTIVE", "PREVIOUSLY_ACTIVE")

.isclose <- function(a, b, rel_tol = 1e-9, abs_tol = 0) {
  if (a == b) return(TRUE)
  if (is.infinite(a) || is.infinite(b)) return(FALSE)
  abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
}

# Matrix exponential for rate matrices with non-negative off-diagonals and
# non-positive row sums (Melloy & Bennett, 1993), as in msprime.
.matrix_exponential <- function(A, n = 5) {
  dA <- -diag(A)
  dmax <- max(dA)
  N <- nrow(A)
  expA <- diag(N)
  if (dmax > 0) {
    P <- A / dmax
    diag(P) <- 1 - dA / dmax
    nscales <- as.integer(max(0, ceiling(log2(dmax) - log2(0.2))))
    t <- dmax / (2^nscales)
    Pk <- diag(N)
    for (k in seq_len(n)) {
      Pk <- (t / k) * (P %*% Pk)
      expA <- expA + Pk
    }
    expA <- expA * exp(-t)
    for (k in seq_len(nscales)) expA <- expA %*% expA
  }
  expA
}

#' Demography debugger
#'
#' Computes and displays the state of the populations in each epoch of a
#' [Demography]. Obtain one with `demography$debug()`.
#'
#' Each epoch is a list with elements `index` (zero-based), `start_time`,
#' `end_time`, `populations` (a list of population states with `id`, `name`,
#' `start_size`, `end_size`, `growth_rate`, `state` and `active`),
#' `migration_matrix` and `events`.
#'
#' @field demography The (resolved) demography.
#' @field num_populations The number of populations.
#' @field epochs The list of epochs.
#' @field num_epochs The number of epochs.
#' @field epoch_start_time The start times of the epochs.
#' @field population_size_history A (num_populations x num_epochs) matrix of
#'   the starting size of each population in each epoch.
#' @export
DemographyDebugger <- R6::R6Class("DemographyDebugger",
  cloneable = FALSE,
  public = list(
    demography = NULL,
    num_populations = NULL,
    epochs = NULL,
    #' @description Create a debugger.
    #' @param Ne,population_configurations,migration_matrix,demographic_events,model
    #'   Deprecated 0.x style model specification.
    #' @param demography A [Demography].
    initialize = function(Ne = 1, population_configurations = NULL,
                          migration_matrix = NULL, demographic_events = NULL,
                          model = NULL, demography = NULL) {
      if (is.null(demography)) {
        demography <- Demography$from_old_style(
          population_configurations, migration_matrix = migration_matrix,
          demographic_events = demographic_events, Ne = Ne, ignore_sample_size = TRUE)
      }
      self$demography <- demography$validate()
      self$num_populations <- demography$num_populations
      private$make_epochs()
      private$check_misspecification()
    },
    #' @description Population sizes at the given times (rows) for each
    #'   population (columns).
    #' @param steps Times ago.
    population_size_trajectory = function(steps) {
      N_t <- matrix(0, length(steps), self$num_populations)
      for (j in seq_along(steps)) N_t[j, ] <- private$pop_size_and_migration_at_t(steps[j])$N
      N_t
    },
    #' @description Lineage location probabilities: an array `P[j, a, b]`
    #'   giving the probability that a lineage sampled in population `a` (at
    #'   `sample_time`) is in population `b` at time `steps[j]` (populations
    #'   are indexed by ID + 1).
    #' @param steps Times at which to compute probabilities.
    #' @param sample_time The sampling time.
    lineage_probabilities = function(steps, sample_time = 0) {
      num_pops <- self$num_populations
      P <- diag(num_pops)
      steps <- as.double(steps)
      epoch_breaks <- setdiff(self$epoch_start_time, steps)
      all_steps <- c(steps, epoch_breaks)
      sampling <- if (!sample_time %in% all_steps) sample_time else numeric(0)
      all_steps <- c(all_steps, sampling)
      ix <- order(all_steps, method = "radix")
      all_steps <- all_steps[ix]
      keep_steps <- c(rep(TRUE, length(steps)), rep(FALSE, length(epoch_breaks)),
                      rep(FALSE, length(sampling)))[ix]
      if (anyDuplicated(all_steps)) stop("steps must be unique", call. = FALSE)
      P_out <- array(0, c(length(all_steps), num_pops, num_pops))
      first_step <- 1
      while (all_steps[first_step] < sample_time) first_step <- first_step + 1
      P_out[first_step, , ] <- P
      mm <- private$lineage_movements()
      if (first_step < length(all_steps)) {
        for (jj in first_step:(length(all_steps) - 1)) {
          t_j <- all_steps[jj]
          if (t_j > sample_time) {
            for (m in mm) {
              if (m$time == t_j) {
                S <- diag(num_pops)
                S[m$source + 1, m$dest + 1] <- m$proportion
                S[m$source + 1, m$source + 1] <- 1 - m$proportion
                P <- P %*% S
              }
            }
          }
          M <- private$pop_size_and_migration_at_t(t_j)$M
          dt <- all_steps[jj + 1] - all_steps[jj]
          dM <- diag(rowSums(M), num_pops)
          P <- P %*% .matrix_exponential(dt * (M - dM))
          P_out[jj + 1, , ] <- P
        }
      }
      P_out[keep_steps, , , drop = FALSE]
    },
    #' @description Epoch intervals in which lineages could be present in each
    #'   population.
    #' @param samples A list of population identifiers or [SampleSet] objects
    #'   (default: all populations with default sampling time 0).
    #' @return A list with `intervals` (a two column matrix of epoch start/end)
    #'   and `indicators` (a logical matrix, one row per interval and one column
    #'   per population).
    possible_lineage_locations = function(samples = NULL) {
      demog <- self$demography
      if (is.null(samples)) {
        samples <- lapply(Filter(function(p) p$default_sampling_time == 0,
                                 demog$populations), function(p) p$id)
      }
      sampling_times <- list()
      for (sample in samples) {
        if (inherits(sample, "SampleSet") || inherits(sample, "Sample")) {
          pop_id <- demog$get(sample$population)$id
          sample_time <- sample$time %||% demog$get(pop_id)$default_sampling_time
          if (inherits(sample, "SampleSet") && sample$num_samples <= 0) next
        } else {
          pop <- demog$get(sample)
          pop_id <- pop$id
          sample_time <- pop$default_sampling_time
        }
        key <- format(sample_time, digits = 17)
        sampling_times[[key]] <- list(time = sample_time,
                                      pops = unique(c(sampling_times[[key]]$pops, pop_id)))
      }
      times <- vapply(sampling_times, function(s) s$time, numeric(1))
      all_steps <- sort(unique(c(self$epoch_start_time, times)))
      epochs <- cbind(all_steps, c(all_steps[-1], Inf))
      all_steps <- c(all_steps, all_steps[length(all_steps)] + 1)
      indicators <- matrix(FALSE, nrow(epochs), self$num_populations)
      for (st in sampling_times) {
        P_out <- self$lineage_probabilities(all_steps, sample_time = st$time)
        for (e in seq_len(nrow(epochs))) {
          if (epochs[e, 2] <= st$time) next
          P <- P_out[e + 1, , , drop = FALSE]
          for (pop_id in st$pops) {
            indicators[e, P[1, pop_id + 1, ] > 0] <- TRUE
          }
        }
      }
      # Join adjacent epochs with identical indicators.
      keep_int <- list()
      keep_ind <- list()
      e <- 1
      while (e <= nrow(epochs)) {
        this <- epochs[e, ]
        f <- e
        while (f + 1 <= nrow(epochs) && all(indicators[f + 1, ] == indicators[e, ])) {
          f <- f + 1
          this[2] <- epochs[f, 2]
        }
        keep_int[[length(keep_int) + 1]] <- this
        keep_ind[[length(keep_ind) + 1]] <- indicators[e, ]
        e <- f + 1
      }
      intervals <- do.call(rbind, keep_int)
      colnames(intervals) <- c("start", "end")
      ind <- do.call(rbind, keep_ind)
      colnames(ind) <- names(demog)
      list(intervals = intervals, indicators = ind)
    },
    #' @description Mean pairwise coalescence time of the specified lineages.
    #' @param lineages Named vector/list mapping populations (names or IDs) to
    #'   numbers of monoploid sample lineages.
    #' @param min_pop_size Smallest allowed population size.
    #' @param steps Initial time discretisation.
    #' @param rtol Relative tolerance.
    #' @param max_iter Maximum number of refinements.
    mean_coalescence_time = function(lineages, min_pop_size = 1, steps = NULL,
                                     rtol = 0.005, max_iter = 12) {
      mean_time <- function(steps, P) {
        dt <- diff(steps)
        dP <- diff(P)
        dlogP <- suppressWarnings(diff(log(P)))
        nz <- dP < 0 & (P[-1] * P[-length(P)]) > 0
        const <- dP == 0
        sum(dt[const] * P[-length(P)][const]) + sum(dt[nz] * dP[nz] / dlogP[nz])
      }
      if (is.null(steps)) {
        last_N <- max(self$population_size_history[, self$num_epochs])
        last_epoch <- self$epoch_start_time[self$num_epochs]
        steps <- sort(unique(c(seq(0, last_epoch + 12 * last_N, length.out = 101),
                               self$epoch_start_time)))
      }
      p_diff <- m_diff <- Inf
      last_P <- Inf
      n <- 0
      while (n < max_iter && (last_P > rtol || p_diff > rtol / 20 || m_diff > rtol / 20)) {
        last_steps <- steps
        if (n == 0) {
          P1 <- self$coalescence_rate_trajectory(last_steps, lineages, min_pop_size,
                                                 double_step_validation = FALSE)$p
          m1 <- mean_time(last_steps, P1)
        }
        if (last_P > rtol) {
          steps <- c(steps, seq(steps[length(steps)], steps[length(steps)] * 1.2,
                                length.out = 20)[-1])
        } else {
          steps <- sort(c(steps, steps[-length(steps)] + diff(steps) / 2))
        }
        P2 <- self$coalescence_rate_trajectory(steps, lineages, min_pop_size,
                                               double_step_validation = FALSE)$p
        m2 <- mean_time(steps, P2)
        keep <- steps %in% last_steps
        p_diff <- max(abs(P1 - P2[keep]))
        m_diff <- abs(m1 - m2) / m2
        last_P <- P2[length(P2)]
        n <- n + 1
        P1 <- P2
        m1 <- m2
      }
      if (n == max_iter) {
        stop("Did not converge on an adequate discretisation: Increase max_iter or rtol.",
             call. = FALSE)
      }
      m2
    },
    #' @description Mean coalescence rates and probabilities of not having
    #'   coalesced, for pairs of the specified lineages at each time in `steps`.
    #' @param steps Increasing, non-negative times.
    #' @param lineages Named vector/list mapping populations to numbers of
    #'   monoploid lineages.
    #' @param min_pop_size Smallest allowed population size.
    #' @param double_step_validation Check the discretisation by doubling it.
    #' @return A list with elements `rates` and `p`.
    coalescence_rate_trajectory = function(steps, lineages, min_pop_size = 1,
                                           double_step_validation = TRUE) {
      num_samples <- numeric(self$num_populations)
      for (key in names(lineages)) {
        pop_id <- .resolve_population_key(key, self$demography)
        num_samples[pop_id + 1] <- num_samples[pop_id + 1] + lineages[[key]]
      }
      steps <- as.double(steps)
      if (length(steps) > 1 && !all(diff(steps) > 0)) {
        stop("`steps` must be a sequence of increasing times.", call. = FALSE)
      }
      if (any(steps < 0)) stop("`steps` must be non-negative", call. = FALSE)
      res <- private$calc_coalescence_rate_trajectory(steps, num_samples, min_pop_size)
      if (double_step_validation && length(steps) > 1) {
        double_steps <- sort(c(steps, steps[-length(steps)] + diff(steps) / 2))
        resd <- private$calc_coalescence_rate_trajectory(double_steps, num_samples,
                                                         min_pop_size)
        idx <- seq(1, length(double_steps), by = 2)
        close <- function(a, b) {
          ok <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) &
                                           abs(a - b) <= 1e-8 + 1e-3 * abs(b))
          all(ok)
        }
        if (!(close(res$rates, resd$rates[idx]) && close(res$p, resd$p[idx]))) {
          warning("Doubling the number of steps has resulted in different predictions, ",
                  "please re-run with smaller step sizes to ensure numerical accuracy.",
                  call. = FALSE)
        }
      }
      res
    },
    #' @description Print the epochs.
    #' @param ... Ignored.
    print = function(...) {
      cat("DemographyDebugger\n")
      for (epoch in self$epochs) {
        if (epoch$index > 0) {
          cat(sprintf("Events @ generation %s\n", .g(epoch$start_time, 3)))
          for (e in epoch$events) {
            cat(sprintf("  - %s: %s\n", .event_type_str[[class(e)[1]]], .event_effect(e)))
          }
        }
        cat(sprintf("Epoch[%d]: [%s, %s) generations\n", epoch$index,
                    .g(epoch$start_time, 3), .g(epoch$end_time, 3)))
        active <- Filter(function(p) p$active, epoch$populations)
        cat(sprintf("  Populations (total=%d active=%d)\n", length(epoch$populations),
                    length(active)))
        if (length(active) > 0) {
          df <- data.frame(
            name = vapply(active, function(p) p$name, character(1)),
            start = vapply(active, function(p) sprintf("%.1f", p$start_size), character(1)),
            end = vapply(active, function(p) sprintf("%.1f", p$end_size), character(1)),
            growth_rate = vapply(active, function(p) .g(p$growth_rate, 3), character(1)),
            stringsAsFactors = FALSE)
          if (length(active) > 1) {
            for (other in active) {
              df[[other$name]] <- vapply(active, function(p) {
                .g(epoch$migration_matrix[p$id + 1, other$id + 1], 4)
              }, character(1))
            }
          }
          print(df, row.names = FALSE)
        }
      }
      invisible(self)
    }
  ),
  active = list(
    num_epochs = function() length(self$epochs),
    epoch_start_time = function() vapply(self$epochs, function(e) e$start_time, numeric(1)),
    population_size_history = function() {
      m <- matrix(0, self$num_populations, self$num_epochs)
      for (j in seq_len(self$num_epochs)) {
        for (k in seq_len(self$num_populations)) {
          m[k, j] <- self$epochs[[j]]$populations[[k]]$start_size
        }
      }
      m
    }
  ),
  private = list(
    make_epochs = function() {
      self$epochs <- list()
      simulator <- .parse_sim_ancestry(demography = self$demography, init_for_debugger = TRUE)
      on.exit(simulator$free(), add = TRUE)
      start_time <- 0
      end_time <- 0
      event_index <- 1
      all_events <- self$demography$events
      while (!is.infinite(end_time)) {
        events <- list()
        while (event_index <= length(all_events) &&
               .isclose(all_events[[event_index]]$time, start_time, abs_tol = 1e-9)) {
          events[[length(events) + 1]] <- all_events[[event_index]]
          event_index <- event_index + 1
        }
        end_time <- simulator$debug_demography()
        state <- simulator$state()
        populations <- lapply(seq_len(self$num_populations), function(j) {
          pc <- state$population_configuration[[j]]
          list(id = j - 1, name = self$demography$populations[[j]]$name,
               start_size = simulator$compute_population_size(j - 1, start_time),
               end_size = simulator$compute_population_size(j - 1, end_time),
               growth_rate = pc$growth_rate, state = .POP_STATE[pc$state + 1],
               active = pc$state == 1)
        })
        self$epochs[[length(self$epochs) + 1]] <- list(
          index = length(self$epochs), start_time = start_time, end_time = end_time,
          populations = populations, migration_matrix = state$migration_matrix,
          events = events)
        start_time <- end_time
      }
    },
    check_misspecification = function() {
      merged <- integer(0)
      for (epoch in self$epochs) {
        for (de in epoch$events) {
          if (inherits(de, "MassMigration") && de$proportion == 1) {
            merged <- unique(c(merged, .convert_id(de$source, self$demography)))
          }
        }
        mm <- epoch$migration_matrix
        for (k in merged) {
          if (any(mm[k + 1, ] != 0) || any(mm[, k + 1] != 0)) {
            .warnf(paste0("Non-zero migration rates exist after merging population %d. ",
                          "This almost certainly indicates demographic misspecification."),
                   k)
          }
        }
      }
    },
    lineage_movements = function() {
      out <- list()
      for (e in self$demography$events) {
        if (inherits(e, "LineageMovementEvent")) {
          for (lm in .as_lineage_movements(e, self$demography)) {
            lm$time <- e$time
            out[[length(out) + 1]] <- lm
          }
        }
      }
      out
    },
    pop_size_and_migration_at_t = function(t) {
      j <- 1
      while (self$epochs[[j]]$end_time <= t) j <- j + 1
      epoch <- self$epochs[[j]]
      N <- vapply(epoch$populations, function(p) {
        p$start_size * exp(-p$growth_rate * (t - epoch$start_time))
      }, numeric(1))
      list(N = N, M = epoch$migration_matrix)
    },
    calc_coalescence_rate_trajectory = function(steps, num_samples, min_pop_size) {
      num_pops <- self$num_populations
      n2 <- num_pops^2
      # IA[x, y] (1-based) indexes the pair state (x, y) in row-major order.
      IA <- matrix(seq_len(n2), num_pops, num_pops, byrow = TRUE)
      P <- matrix(0, 1, n2)
      for (x in seq_len(num_pops)) {
        for (y in seq_len(num_pops)) {
          P[1, IA[x, y]] <- num_samples[x] * (num_samples[y] - (x == y))
        }
      }
      P <- P / sum(P)
      # Unlike msprime, time 0 is only added once (it may already be in steps).
      epoch_breaks <- unique(setdiff(c(0, self$epoch_start_time), steps))
      steps_b <- c(steps, epoch_breaks)
      ix <- order(steps_b, method = "radix")
      steps_b <- steps_b[ix]
      keep_steps <- c(rep(TRUE, length(steps)), rep(FALSE, length(epoch_breaks)))[ix]
      mm <- private$lineage_movements()
      num_steps <- length(steps_b)
      Identity <- diag(num_pops)
      r <- numeric(num_steps)
      p_t <- numeric(num_steps)
      for (j in seq_len(num_steps)) {
        time <- steps_b[j]
        NM <- private$pop_size_and_migration_at_t(time)
        N <- NM$N
        M <- NM$M
        C <- matrix(0, n2, n2)
        for (idx in seq_len(num_pops)) {
          C[IA[idx, idx], IA[idx, idx]] <- 1 / (2 * max(min_pop_size, N[idx]))
        }
        dM <- diag(rowSums(M), num_pops)
        for (m in mm) {
          if (m$time == time) {
            a <- m$source + 1
            b <- m$dest + 1
            p <- m$proportion
            S <- diag(n2)
            for (x in seq_len(num_pops)) {
              if (x == a) {
                S[IA[a, a], IA[a, b]] <- p * (1 - p)
                S[IA[a, a], IA[b, a]] <- p * (1 - p)
                S[IA[a, a], IA[b, b]] <- p^2
                S[IA[a, a], IA[a, a]] <- (1 - p)^2
              } else {
                S[IA[x, a], IA[x, b]] <- p
                S[IA[a, x], IA[b, x]] <- p
                S[IA[x, a], IA[x, a]] <- 1 - p
                S[IA[a, x], IA[a, x]] <- 1 - p
              }
            }
            P <- P %*% S
          }
        }
        p_notcoal <- sum(P)
        p_t[j] <- p_notcoal
        r[j] <- if (p_notcoal > 0) sum(P %*% C) / p_notcoal else NaN
        if (j < num_steps) {
          dt <- steps_b[j + 1] - steps_b[j]
          G <- (kronecker(M - dM, Identity) + kronecker(Identity, M - dM)) - C
          P <- P %*% .matrix_exponential(dt * G)
        }
      }
      list(rates = r[keep_steps], p = p_t[keep_steps])
    }
  )
)

.resolve_population_key <- function(key, demography) {
  if (key %in% names(demography)) return(demography$get(key)$id)
  id <- suppressWarnings(as.numeric(key))
  if (!is.na(id)) return(demography$get(id)$id)
  .stopf("Population with name '%s' not found", key)
}

.assert_equivalent <- function(self, other, rel_tol, abs_tol) {
  fail <- function(...) {
    stop(structure(class = c("msprimeR_not_equivalent", "error", "condition"),
                   list(message = sprintf(...), call = NULL)))
  }
  if (!inherits(other, "Demography")) stop("Not a Demography", call. = FALSE)
  if (self$num_populations != other$num_populations) {
    fail("Number of populations not equal: %d != %d", self$num_populations,
         other$num_populations)
  }
  self_dbg <- self$debug()
  other_dbg <- other$debug()
  for (j in seq_len(self$num_populations)) {
    sp <- self_dbg$demography$populations[[j]]
    op <- other_dbg$demography$populations[[j]]
    if (sp$name != op$name) fail("Population names differ: %s != %s", sp$name, op$name)
    if (!.isclose(sp$default_sampling_time %||% 0, op$default_sampling_time %||% 0)) {
      fail("Sampling times not equal for %s", sp$name)
    }
  }
  if (self_dbg$num_epochs != other_dbg$num_epochs) {
    fail("Number of epochs not equal: %d != %d", self_dbg$num_epochs, other_dbg$num_epochs)
  }
  for (j in seq_len(self_dbg$num_epochs)) {
    se <- self_dbg$epochs[[j]]
    oe <- other_dbg$epochs[[j]]
    if (!.isclose(se$start_time, oe$start_time, rel_tol, abs_tol)) {
      fail("Epoch[%d] at different times: %s != %s", j - 1, se$start_time, oe$start_time)
    }
    for (k in seq_along(se$populations)) {
      sp <- se$populations[[k]]
      op <- oe$populations[[k]]
      if (sp$state != op$state) {
        fail("State mismatch in populations in epoch[%d], %s", j - 1, sp$name)
      }
      if (sp$active) {
        if (!.isclose(sp$start_size, op$start_size, rel_tol, abs_tol)) {
          fail("Population start_size not equal in epoch[%d], %s: %s != %s", j - 1,
               sp$name, sp$start_size, op$start_size)
        }
        if (!.isclose(sp$growth_rate, op$growth_rate, rel_tol, abs_tol)) {
          fail("Population growth_rate not equal in epoch[%d], %s: %s != %s", j - 1,
               sp$name, sp$growth_rate, op$growth_rate)
        }
      }
    }
    a <- se$migration_matrix
    b <- oe$migration_matrix
    if (!all(abs(a - b) <= abs_tol + rel_tol * abs(b))) {
      fail("Migration matrices in epoch[%d] not equal", j - 1)
    }
    s_lm <- .normalise_lineage_movements(se$events, self)
    o_lm <- .normalise_lineage_movements(oe$events, other)
    if (!setequal(names(s_lm), names(o_lm))) {
      fail("Mismatch in the set of populations affected by lineage movements")
    }
    for (src in names(s_lm)) {
      a <- s_lm[[src]]
      b <- o_lm[[src]]
      if (length(a) != length(b)) {
        fail("Mismatch in number of normalised lineage movements out of %s", src)
      }
      for (q in seq_along(a)) {
        if (a[[q]]$dest != b[[q]]$dest) fail("Mismatch in lineage movement destination")
        if (!.isclose(a[[q]]$proportion, b[[q]]$proportion, rel_tol, abs_tol)) {
          fail("Mismatch in normalised lineage movement proportions from %s", src)
        }
      }
    }
    if (any(vapply(c(se$events, oe$events), inherits, logical(1), "StateChangeEvent"))) {
      stop("State change events not currently supported in equivalent.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

.normalise_lineage_movements <- function(events, demography) {
  ret <- list()
  for (event in events) {
    if (inherits(event, "LineageMovementEvent")) {
      for (lm in .as_lineage_movements(event, demography)) {
        key <- as.character(lm$source)
        ret[[key]] <- c(ret[[key]], list(lm))
      }
    }
  }
  for (key in names(ret)) {
    P <- .sequential_to_proportions(vapply(ret[[key]], function(m) m$proportion, numeric(1)))
    dests <- vapply(ret[[key]], function(m) m$dest, numeric(1))
    o <- order(dests, P)
    S <- .proportions_to_sequential(P[o])
    ret[[key]] <- lapply(seq_along(o), function(j) {
      list(source = as.numeric(key), dest = dests[o][j], proportion = S[j])
    })
  }
  ret
}
