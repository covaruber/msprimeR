# Mutation models and sim_mutations() (port of msprime/mutations.py).

.ACGT_ALLELES <- c("A", "C", "G", "T")
.AMINO_ACIDS <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F",
                  "P", "S", "T", "W", "Y", "V")

#' @rdname mutation_models
#' @export
BINARY <- 0
#' @rdname mutation_models
#' @export
NUCLEOTIDES <- 1

MutationModelClass <- R6::R6Class("MutationModel",
  public = list(
    args = NULL,
    asdict = function() self$args,
    ll_spec = function() stop("Abstract mutation model", call. = FALSE),
    update_state = function(state) invisible(NULL)
  )
)

MatrixMutationModelClass <- R6::R6Class("MatrixMutationModel",
  inherit = MutationModelClass,
  public = list(
    alleles = NULL,
    root_distribution = NULL,
    transition_matrix = NULL,
    initialize = function(alleles, root_distribution, transition_matrix, args = NULL) {
      alleles <- as.character(unlist(alleles))
      root_distribution <- as.double(unlist(root_distribution))
      if (is.list(transition_matrix)) transition_matrix <- do.call(rbind, transition_matrix)
      transition_matrix <- as.matrix(transition_matrix)
      storage.mode(transition_matrix) <- "double"
      .Call(msprimer_check_matrix_model, alleles, root_distribution, transition_matrix)
      self$alleles <- alleles
      self$root_distribution <- root_distribution
      self$transition_matrix <- unname(transition_matrix)
      self$args <- args %||% list(alleles = alleles, root_distribution = root_distribution,
                                  transition_matrix = transition_matrix)
    },
    ll_spec = function() {
      list(type = "matrix", alleles = self$alleles,
           root_distribution = self$root_distribution,
           transition_matrix = self$transition_matrix)
    },
    print = function(...) {
      cat(sprintf("Mutation model with alleles [%s]\n",
                  paste0("'", self$alleles, "'", collapse = ", ")))
      cat("  root distribution:", paste(format(self$root_distribution), collapse = " "), "\n")
      cat("  transition matrix:\n")
      for (j in seq_len(nrow(self$transition_matrix))) {
        cat("    ", paste(format(self$transition_matrix[j, ]), collapse = " "), "\n")
      }
      invisible(self)
    }
  )
)

.matrix_model_subclass <- function(name) {
  R6::R6Class(name, inherit = MatrixMutationModelClass)
}

SLiMv5MutationModelClass <- R6::R6Class("SLiMv5MutationModel",
  inherit = MutationModelClass,
  public = list(
    type = NULL, next_id = NULL, slim_generation = NULL, block_size = NULL,
    initialize = function(type, next_id = 0, slim_generation = 1, block_size = 0) {
      if (!.isinteger(type)) stop("type must be an integer", call. = FALSE)
      if (type < 0) stop("type must be >= 0", call. = FALSE)
      if (!.isinteger(next_id) || next_id < 0) stop("next_id must be >= 0", call. = FALSE)
      self$type <- type
      self$next_id <- next_id
      self$slim_generation <- slim_generation
      self$block_size <- block_size
      self$args <- list(type = type, next_id = next_id, slim_generation = slim_generation,
                        block_size = block_size)
    },
    ll_spec = function() {
      list(type = "slim_v5", type_id = self$type, next_id = self$next_id,
           slim_generation = self$slim_generation, block_size = self$block_size)
    },
    update_state = function(state) {
      self$next_id <- state$next_id
      self$args$next_id <- state$next_id
    },
    print = function(...) {
      cat(sprintf("Mutation model for SLiM mutations of type m%s\n  next ID: %s\n",
                  format(self$type), format(self$next_id)))
      invisible(self)
    }
  )
)

SLiMv6MutationModelClass <- R6::R6Class("SLiMv6MutationModel",
  inherit = MutationModelClass,
  public = list(
    next_id = NULL, block_size = NULL,
    initialize = function(next_id = 0, block_size = 0) {
      if (!.isinteger(next_id) || next_id < 0) stop("next_id must be >= 0", call. = FALSE)
      self$next_id <- next_id
      self$block_size <- block_size
      self$args <- list(next_id = next_id, block_size = block_size)
    },
    ll_spec = function() {
      list(type = "slim_v6", next_id = self$next_id, block_size = self$block_size)
    },
    update_state = function(state) {
      self$next_id <- state$next_id
      self$args$next_id <- state$next_id
    },
    print = function(...) {
      cat(sprintf("Mutation model for SLiM mutations, v6+. Next ID: %s\n",
                  format(self$next_id)))
      invisible(self)
    }
  )
)

InfiniteAllelesClass <- R6::R6Class("InfiniteAlleles",
  inherit = MutationModelClass,
  public = list(
    start_allele = NULL, next_allele = NULL,
    initialize = function(start_allele = 0) {
      if (!.isinteger(start_allele) || start_allele < 0) {
        stop("start_allele must be a non-negative integer", call. = FALSE)
      }
      self$start_allele <- start_allele
      self$next_allele <- start_allele
      self$args <- list(start_allele = start_allele)
    },
    ll_spec = function() {
      list(type = "infinite_alleles", start_allele = self$start_allele,
           next_allele = self$next_allele)
    },
    update_state = function(state) self$next_allele <- state$next_allele,
    print = function(...) {
      cat(sprintf(paste0("Infinite alleles mutation model, beginning with allele %s\n",
                         "    next allele: %s\n"), format(self$start_allele),
                  format(self$next_allele)))
      invisible(self)
    }
  )
)

#' Mutation models
#'
#' Models of sequence evolution for [sim_mutations()]. Matrix mutation models
#' have a finite set of `alleles`, a `root_distribution` from which ancestral
#' states are drawn and a `transition_matrix` giving the probabilities of each
#' derived state given the parent state.
#'
#' * `MatrixMutationModel(alleles, root_distribution, transition_matrix)`: a
#'   general matrix model.
#' * `BinaryMutationModel()`: alleles "0" and "1".
#' * `JC69()`, `HKY(kappa)`, `F84(kappa)`, `GTR(relative_rates)`: nucleotide
#'   models.
#' * `BLOSUM62()`, `PAM()`: amino acid models.
#' * `MicrosatMutationModel()`, `SMM()`, `TPM()`, `EL2()`: microsatellite
#'   models (Sainudiin et al. 2004).
#' * `InfiniteAlleles(start_allele)`: each new allele is a new integer.
#' * `SLiMv5MutationModel(type)`, `SLiMv6MutationModel()`: SLiM-compatible
#'   mutations; `SLiMMutationModel()` is a deprecated alias of the v5 model.
#' * `InfiniteSites(alphabet)`: legacy model for [mutate()].
#'
#' Stateful models (`InfiniteAlleles`, SLiM models) are updated in place by
#' [sim_mutations()], so that subsequent calls continue the allele numbering,
#' as in Python's msprime.
#'
#' @param alleles A character vector of alleles.
#' @param root_distribution Probabilities of each allele at the root.
#' @param transition_matrix A square matrix of transition probabilities (rows
#'   sum to one).
#' @param state_independent If TRUE, mutations are generated independently of
#'   the previous state (including silent mutations).
#' @param kappa Transition/transversion rate ratio.
#' @param equilibrium_frequencies Equilibrium allele frequencies (A, C, G, T).
#' @param relative_rates GTR rates in the order A<>C, A<>G, A<>T, C<>G, C<>T,
#'   G<>T.
#' @param s,u,v,p,m Microsatellite model parameters.
#' @param lo,hi Lower and upper bounds on the repeat number.
#' @param start_allele The first allele for the infinite alleles model.
#' @param type The SLiM mutation type.
#' @param next_id The next SLiM mutation ID.
#' @param slim_generation The SLiM generation.
#' @param block_size Allocation block size for derived states.
#' @param alphabet `BINARY` or `NUCLEOTIDES`.
#' @return A mutation model object.
#' @name mutation_models
NULL

#' @rdname mutation_models
#' @export
MatrixMutationModel <- function(alleles, root_distribution, transition_matrix) {
  MatrixMutationModelClass$new(alleles, root_distribution, transition_matrix)
}

.BinaryClass <- .matrix_model_subclass("BinaryMutationModel")
.JC69Class <- .matrix_model_subclass("JC69")
.HKYClass <- .matrix_model_subclass("HKY")
.F84Class <- .matrix_model_subclass("F84")
.GTRClass <- .matrix_model_subclass("GTR")
.BLOSUM62Class <- .matrix_model_subclass("BLOSUM62")
.PAMClass <- .matrix_model_subclass("PAM")
.MicrosatClass <- .matrix_model_subclass("MicrosatMutationModel")
.SMMClass <- R6::R6Class("SMM", inherit = .MicrosatClass)
.TPMClass <- R6::R6Class("TPM", inherit = .MicrosatClass)
.EL2Class <- R6::R6Class("EL2", inherit = .MicrosatClass)
.InfiniteSitesClass <- .matrix_model_subclass("InfiniteSites")

#' @rdname mutation_models
#' @export
BinaryMutationModel <- function(state_independent = FALSE) {
  tm <- if (state_independent) matrix(0.5, 2, 2) else matrix(c(0, 1, 1, 0), 2, 2)
  .BinaryClass$new(c("0", "1"), c(1, 0), tm,
                   args = list(state_independent = state_independent))
}

#' @rdname mutation_models
#' @export
JC69 <- function(state_independent = FALSE) {
  if (state_independent) {
    tm <- matrix(1 / 4, 4, 4)
  } else {
    tm <- matrix(1 / 3, 4, 4)
    diag(tm) <- 0
  }
  .JC69Class$new(.ACGT_ALLELES, rep(0.25, 4), tm,
                 args = list(state_independent = state_independent))
}

.scale_columns <- function(m, v) m * rep(v, each = nrow(m))

#' @rdname mutation_models
#' @export
HKY <- function(kappa, equilibrium_frequencies = NULL, root_distribution = NULL,
                state_independent = FALSE) {
  args <- list(kappa = kappa, equilibrium_frequencies = equilibrium_frequencies,
               root_distribution = root_distribution, state_independent = state_independent)
  pi <- if (is.null(equilibrium_frequencies)) rep(0.25, 4) else as.double(equilibrium_frequencies)
  root <- root_distribution %||% pi
  tm <- matrix(1.0, 4, 4)
  tm[cbind(c(1, 2, 3, 4), c(3, 4, 1, 2))] <- kappa
  tm <- .scale_columns(tm, pi)
  diag(tm) <- if (state_independent) kappa * pi else 0
  tm <- tm / max(.np_row_sums(tm))
  diag(tm) <- 1.0 - (.np_row_sums(tm) - diag(tm))
  .HKYClass$new(.ACGT_ALLELES, root, tm, args = args)
}

#' @rdname mutation_models
#' @export
F84 <- function(kappa, equilibrium_frequencies = NULL, root_distribution = NULL,
                state_independent = FALSE) {
  args <- list(kappa = kappa, equilibrium_frequencies = equilibrium_frequencies,
               root_distribution = root_distribution, state_independent = state_independent)
  pi <- if (is.null(equilibrium_frequencies)) rep(0.25, 4) else as.double(equilibrium_frequencies)
  root <- root_distribution %||% pi
  tm <- matrix(1.0, 4, 4)
  p_AG <- pi[1] + pi[3]
  p_CT <- pi[2] + pi[4]
  gamma <- 1 + (kappa - 1) / c(p_AG, p_CT)
  tm[1, 3] <- tm[3, 1] <- gamma[1]
  tm[2, 4] <- tm[4, 2] <- gamma[2]
  tm <- .scale_columns(tm, pi)
  diag(tm) <- if (state_independent) pi * c(gamma, gamma) else 0
  tm <- tm / max(.np_row_sums(tm))
  diag(tm) <- 1.0 - (.np_row_sums(tm) - diag(tm))
  .F84Class$new(.ACGT_ALLELES, root, tm, args = args)
}

#' @rdname mutation_models
#' @export
GTR <- function(relative_rates, equilibrium_frequencies = NULL, root_distribution = NULL,
                state_independent = FALSE) {
  args <- list(relative_rates = relative_rates,
               equilibrium_frequencies = equilibrium_frequencies,
               root_distribution = root_distribution, state_independent = state_independent)
  if (length(relative_rates) != 6) stop("relative_rates must have 6 entries", call. = FALSE)
  pi <- if (is.null(equilibrium_frequencies)) rep(0.25, 4) else as.double(equilibrium_frequencies)
  root <- root_distribution %||% pi
  tm <- matrix(0, 4, 4)
  upper <- cbind(c(1, 1, 1, 2, 2, 3), c(2, 3, 4, 3, 4, 4))
  tm[upper] <- as.double(relative_rates)
  tm <- tm + t(tm)
  if (state_independent) diag(tm) <- 1.0
  tm <- .scale_columns(tm, pi)
  tm <- tm / max(.np_row_sums(tm))
  diag(tm) <- diag(tm) + 1.0 - .np_row_sums(tm)
  .GTRClass$new(.ACGT_ALLELES, root, tm, args = args)
}

.amino_acid_matrix <- function(relative_rates, root) {
  n <- 20
  tm <- matrix(0, n, n)
  idx <- do.call(rbind, lapply(2:n, function(i) cbind(i, seq_len(i - 1))))
  tm[idx] <- relative_rates
  tm <- tm + t(tm)
  tm <- .scale_columns(tm, root)
  tm <- tm / max(.np_row_sums(tm))
  diag(tm) <- 1.0 - .np_row_sums(tm)
  tm
}

#' @rdname mutation_models
#' @export
BLOSUM62 <- function() {
  root <- c(0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074, 0.026, 0.068, 0.099,
            0.058, 0.025, 0.047, 0.039, 0.057, 0.051, 0.013, 0.032, 0.073)
  rr <- c(
    0.735790389698, 0.485391055466, 1.297446705134, 0.543161820899, 0.500964408555,
    3.180100048216, 1.459995310470, 0.227826574209, 0.397358949897, 0.240836614802,
    1.199705704602, 3.020833610064, 1.839216146992, 1.190945703396, 0.329801504630,
    1.170949042800, 1.360574190420, 1.240488508640, 3.761625208368, 0.140748891814,
    5.528919177928, 1.955883574960, 0.418763308518, 1.355872344485, 0.798473248968,
    0.418203192284, 0.609846305383, 0.423579992176, 0.716241444998, 1.456141166336,
    2.414501434208, 0.778142664022, 0.354058109831, 2.435341131140, 1.626891056982,
    0.539859124954, 0.605899003687, 0.232036445142, 0.283017326278, 0.418555732462,
    0.774894022794, 0.236202451204, 0.186848046932, 0.189296292376, 0.252718447885,
    0.800016530518, 0.622711669692, 0.211888159615, 0.218131577594, 0.831842640142,
    0.580737093181, 0.372625175087, 0.217721159236, 0.348072209797, 3.890963773304,
    1.295201266783, 5.411115141489, 1.593137043457, 1.032447924952, 0.285078800906,
    3.945277674515, 2.802427151679, 0.752042440303, 1.022507035889, 0.406193586642,
    0.445570274261, 1.253758266664, 0.983692987457, 0.648441278787, 0.222621897958,
    0.767688823480, 2.494896077113, 0.555415397470, 0.459436173579, 0.984311525359,
    3.364797763104, 6.030559379572, 1.073061184332, 0.492964679748, 0.371644693209,
    0.354861249223, 0.281730694207, 0.441337471187, 0.144356959750, 0.291409084165,
    0.368166464453, 0.714533703928, 1.517359325954, 2.064839703237, 0.266924750511,
    1.773855168830, 1.173275900924, 0.448133661718, 0.494887043702, 0.730628272998,
    0.356008498769, 0.858570575674, 0.926563934846, 0.504086599527, 0.527007339151,
    0.388355409206, 0.374555687471, 1.047383450722, 0.454123625103, 0.233597909629,
    4.325092687057, 1.122783104210, 2.904101656456, 1.582754142065, 1.197188415094,
    1.934870924596, 1.769893238937, 1.509326253224, 1.117029762910, 0.357544412460,
    0.352969184527, 1.752165917819, 0.918723415746, 0.540027644824, 1.169129577716,
    1.729178019485, 0.914665954563, 1.898173634533, 0.934187509431, 1.119831358516,
    1.277480294596, 1.071097236007, 0.641436011405, 0.585407090225, 1.179091197260,
    0.915259857694, 1.303875200799, 1.488548053722, 0.488206118793, 1.005451683149,
    5.151556292270, 0.465839367725, 0.426382310122, 0.191482046247, 0.145345046279,
    0.527664418872, 0.758653808642, 0.407635648938, 0.508358924638, 0.301248600780,
    0.341985787540, 0.691474634600, 0.332243040634, 0.888101098152, 2.074324893497,
    0.252214830027, 0.387925622098, 0.513128126891, 0.718206697586, 0.720517441216,
    0.538222519037, 0.261422208965, 0.470237733696, 0.958989742850, 0.596719300346,
    0.308055737035, 4.218953969389, 0.674617093228, 0.811245856323, 0.717993486900,
    0.951682162246, 6.747260430801, 0.369405319355, 0.796751520761, 0.801010243199,
    4.054419006558, 2.187774522005, 0.438388343772, 0.312858797993, 0.258129289418,
    1.116352478606, 0.530785790125, 0.524253846338, 0.253340790190, 0.201555971750,
    8.311839405458, 2.231405688913, 0.498138475304, 2.575850755315, 0.838119610178,
    0.496908410676, 0.561925457442, 2.253074051176, 0.266508731426, 1.000000000000)
  .BLOSUM62Class$new(.AMINO_ACIDS, root, .amino_acid_matrix(rr, root), args = list())
}

#' @rdname mutation_models
#' @export
PAM <- function() {
  root <- c(0.087127, 0.040904, 0.040432, 0.046872, 0.033474, 0.038255, 0.049530,
            0.088612, 0.033619, 0.036886, 0.085357, 0.080481, 0.014753, 0.039772,
            0.050680, 0.069577, 0.058542, 0.010494, 0.029916, 0.064717)
  rr <- c(
    0.267828, 0.984474, 0.327059, 1.199805, 0.000000, 8.931515, 0.360016, 0.232374,
    0.000000, 0.000000, 0.887753, 2.439939, 1.028509, 1.348551, 0.000000, 1.961167,
    0.000000, 1.493409, 11.388659, 0.000000, 7.086022, 2.386111, 0.087791, 1.385352,
    1.240981, 0.107278, 0.281581, 0.811907, 0.228116, 2.383148, 5.290024, 0.868241,
    0.282729, 6.011613, 0.439469, 0.106802, 0.653416, 0.632629, 0.768024, 0.239248,
    0.438074, 0.180393, 0.609526, 0.000000, 0.076981, 0.406431, 0.154924, 0.341113,
    0.000000, 0.000000, 0.730772, 0.112880, 0.071514, 0.443504, 2.556685, 0.258635,
    4.610124, 3.148371, 0.716913, 0.000000, 1.519078, 0.830078, 0.267683, 0.270475,
    0.460857, 0.180629, 0.717840, 0.896321, 0.000000, 0.000000, 0.000000, 1.127499,
    0.304803, 0.170372, 0.000000, 3.332732, 5.230115, 2.411739, 0.183641, 0.136906,
    0.138503, 0.000000, 0.000000, 0.000000, 0.000000, 0.153478, 0.475927, 1.951951,
    1.565160, 0.000000, 0.921860, 2.485920, 1.028313, 0.419244, 0.133940, 0.187550,
    1.526188, 0.507003, 0.347153, 0.933709, 0.119152, 0.316258, 0.335419, 0.170205,
    0.110506, 4.051870, 1.531590, 4.885892, 0.956097, 1.598356, 0.561828, 0.793999,
    2.322243, 0.353643, 0.247955, 0.171432, 0.954557, 0.619951, 0.459901, 2.427202,
    3.680365, 0.265745, 2.271697, 0.660930, 0.162366, 0.525651, 0.340156, 0.306662,
    0.226333, 1.900739, 0.331090, 1.350599, 1.031534, 0.136655, 0.782857, 5.436674,
    0.000000, 2.001375, 0.224968, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.270564, 0.000000, 0.461776, 0.000000, 0.000000, 0.762354, 0.000000, 0.740819,
    0.000000, 0.244139, 0.078012, 0.946940, 0.000000, 0.953164, 0.000000, 0.214717,
    0.000000, 1.265400, 0.374834, 0.286572, 0.132142, 0.000000, 6.952629, 0.000000,
    0.336289, 0.417839, 0.608070, 2.059564, 0.240368, 0.158067, 0.178316, 0.484678,
    0.346983, 0.367250, 0.538165, 0.438715, 8.810038, 1.745156, 0.103850, 2.565955,
    0.123606, 0.485026, 0.303836, 1.561997, 0.000000, 0.279379)
  .PAMClass$new(.AMINO_ACIDS, root, .amino_acid_matrix(rr, root), args = list())
}

#' Microsatellite rate matrix
#'
#' The jump-process transition matrix of the Sainudiin et al. (2004)
#' microsatellite model, as used by [MicrosatMutationModel()].
#'
#' @param s,u,v,p,m Model parameters (see [mutation_models]).
#' @param lo,hi Bounds on the repeat number.
#' @return A square matrix with `hi - lo + 1` rows.
#' @export
general_microsat_rate_matrix <- function(s, u, v, p, m, lo, hi) {
  alpha_ms <- function(i) max(0, min(1, u - (v * (i - lo))))
  beta_ms <- function(i) 1 + ((i - lo) * s)
  gamma_ms <- function(i, j) {
    num <- m * (1 - m)^(abs(i - j) - 1)
    if (lo <= i && i < j && j <= hi) {
      denom <- 1 - (1 - m)^(hi - i)
    } else if (lo <= j && j < i && i <= hi) {
      denom <- 1 - (1 - m)^(i - lo)
    } else {
      stop("microsat model error. i and j must be between lo and hi", call. = FALSE)
    }
    num / denom
  }
  n <- hi - lo + 1
  rm <- matrix(0, n, n)
  for (i in lo:hi) {
    for (j in lo:hi) {
      ii <- i - lo + 1
      jj <- j - lo + 1
      rm[ii, jj] <- if (i == j - 1) {
        beta_ms(i) * alpha_ms(i) * (p + ((1 - p) * gamma_ms(i, j)))
      } else if (i < j - 1) {
        beta_ms(i) * alpha_ms(i) * (1 - p) * gamma_ms(i, j)
      } else if (i == j + 1) {
        beta_ms(i) * (1 - alpha_ms(i)) * (p + ((1 - p) * gamma_ms(i, j)))
      } else if (i > j + 1) {
        beta_ms(i) * (1 - alpha_ms(i)) * (1 - p) * gamma_ms(i, j)
      } else {
        0
      }
    }
  }
  rm <- rm / max(.np_row_sums(rm))
  diag(rm) <- pmax(0.0, 1.0 - .np_row_sums(rm))
  rm
}

.microsat <- function(cls, s = NULL, u = NULL, v = NULL, p = NULL, m = NULL, lo = NULL,
                      hi = NULL, root_distribution = NULL, args) {
  s <- s %||% 0.0
  u <- u %||% 0.5
  v <- v %||% 0.0
  p <- p %||% 1.0
  m <- m %||% 1.0
  lo <- lo %||% 1
  hi <- hi %||% 50
  params <- list(s, u, v, p, m, lo, hi)
  if (!all(vapply(params, function(x) is.numeric(x) && length(x) == 1, logical(1)))) {
    stop("All parameters must be numeric", call. = FALSE)
  }
  if (!all(vapply(params, is.finite, logical(1)))) {
    stop("All parameters must be finite", call. = FALSE)
  }
  lo <- as.integer(lo)
  hi <- as.integer(hi)
  if (lo > hi) stop("lo > hi", call. = FALSE)
  if (lo < 0) stop("lo must be >= 0", call. = FALSE)
  if (s <= -1 / (hi - lo + 1)) stop("s must be > -1/(hi-lo+1)", call. = FALSE)
  if (u <= 0 || u >= 1) stop("u must be between 0 and 1", call. = FALSE)
  alleles <- as.character(lo:hi)
  tm <- general_microsat_rate_matrix(s, u, v, p, m, lo, hi)
  if (is.null(root_distribution)) {
    e <- eigen(t(tm))
    k <- which(abs(e$values - 1.0) < 1e-8)[1]
    stationary <- Re(e$vectors[, k])
    root_distribution <- stationary / sum(stationary)
  }
  cls$new(alleles, root_distribution, tm, args = args)
}

#' @rdname mutation_models
#' @export
MicrosatMutationModel <- function(s = NULL, u = NULL, v = NULL, p = NULL, m = NULL,
                                  lo = NULL, hi = NULL, root_distribution = NULL) {
  args <- list(s = s, u = u, v = v, p = p, m = m, lo = lo, hi = hi,
               root_distribution = root_distribution)
  .microsat(.MicrosatClass, s, u, v, p, m, lo, hi, root_distribution, args)
}

#' @rdname mutation_models
#' @export
SMM <- function(lo = NULL, hi = NULL, root_distribution = NULL) {
  .microsat(.SMMClass, s = 0, u = 0.5, v = 0, p = 0, m = 1, lo = lo, hi = hi,
            root_distribution = root_distribution,
            args = list(lo = lo, hi = hi, root_distribution = root_distribution))
}

#' @rdname mutation_models
#' @export
TPM <- function(p, m, lo = NULL, hi = NULL, root_distribution = NULL) {
  if (!(0 < p && p < 1)) stop("p must be between 0 and 1", call. = FALSE)
  if (!(0 < m && m < 1)) stop("m must be between 0 and 1", call. = FALSE)
  .microsat(.TPMClass, s = 0, u = 0.5, v = 0, p = p, m = m, lo = lo, hi = hi,
            root_distribution = root_distribution,
            args = list(p = p, m = m, lo = lo, hi = hi, root_distribution = root_distribution))
}

#' @rdname mutation_models
#' @export
EL2 <- function(m, u, v, lo = NULL, hi = NULL, root_distribution = NULL) {
  if (!(0 < m && m < 1)) stop("m must be between 0 and 1", call. = FALSE)
  if (!(0 < u && u < 1)) stop("u must be between 0 and 1", call. = FALSE)
  .microsat(.EL2Class, s = 0, u = u, v = v, p = 0, m = m, lo = lo, hi = hi,
            root_distribution = root_distribution,
            args = list(m = m, u = u, v = v, lo = lo, hi = hi,
                        root_distribution = root_distribution))
}

#' @rdname mutation_models
#' @export
InfiniteSites <- function(alphabet = BINARY) {
  model <- if (identical(alphabet, BINARY)) {
    BinaryMutationModel()
  } else if (identical(alphabet, NUCLEOTIDES)) {
    JC69()
  } else {
    stop("Bad alphabet", call. = FALSE)
  }
  .InfiniteSitesClass$new(model$alleles, model$root_distribution, model$transition_matrix,
                          args = list(alphabet = alphabet))
}

#' @rdname mutation_models
#' @export
InfiniteAlleles <- function(start_allele = 0) InfiniteAllelesClass$new(start_allele)

#' @rdname mutation_models
#' @export
SLiMv5MutationModel <- function(type, next_id = 0, slim_generation = 1, block_size = 0) {
  SLiMv5MutationModelClass$new(type, next_id, slim_generation, block_size)
}

#' @rdname mutation_models
#' @export
SLiMMutationModel <- function(type, next_id = 0, slim_generation = 1, block_size = 0) {
  warning("SLiMMutationModel will be removed in future versions. Use SLiMv5MutationModel ",
          "or SLiMv6MutationModel instead.", call. = FALSE)
  SLiMv5MutationModelClass$new(type, next_id, slim_generation, block_size)
}

#' @rdname mutation_models
#' @export
SLiMv6MutationModel <- function(next_id = 0, block_size = 0) {
  SLiMv6MutationModelClass$new(next_id, block_size)
}

.MUTATION_MODEL_MAP <- list(infinite_alleles = function() InfiniteAlleles(),
                            binary = function() BinaryMutationModel(),
                            jc69 = function() JC69(), blosum62 = function() BLOSUM62(),
                            pam = function() PAM())

.mutation_model_factory <- function(model) {
  if (is.null(model)) return(JC69())
  if (is.character(model) && length(model) == 1) {
    lower <- tolower(model)
    if (!lower %in% names(.MUTATION_MODEL_MAP)) {
      .stopf("Model '%s' unknown. Choose from %s", model,
             paste0("['", paste(sort(names(.MUTATION_MODEL_MAP)), collapse = "', '"), "']"))
    }
    return(.MUTATION_MODEL_MAP[[lower]]())
  }
  if (inherits(model, "MutationModel")) return(model)
  stop("Mutation model must be a string or an instance of MutationModel", call. = FALSE)
}

#' Simulate mutations
#'
#' Simulates mutations on the given tree sequence at the specified rate (per
#' unit of sequence length per unit of time), returning a new tree sequence.
#' The R equivalent of Python's `msprime.sim_mutations()`.
#'
#' @param tree_sequence The tree sequence (see [msprime_ts]).
#' @param rate A single rate or a [RateMap] (default 0).
#' @param random_seed The random seed (see [sim_ancestry()]).
#' @param model A mutation model or model name (default [JC69()]); see
#'   [mutation_models].
#' @param start_time,end_time Restrict mutations to this time interval.
#' @param discrete_genome If TRUE (default), mutations occur at integer
#'   positions and may stack at the same site.
#' @param keep If TRUE (default), keep existing mutations.
#' @param record_provenance If TRUE (default), record a provenance entry.
#' @return A tree sequence.
#' @examples
#' ts <- sim_ancestry(5, sequence_length = 1e4, population_size = 1e4,
#'                    random_seed = 1)
#' mts <- sim_mutations(ts, rate = 1e-7, random_seed = 2)
#' @export
sim_mutations <- function(tree_sequence, rate = NULL, random_seed = NULL, model = NULL,
                          start_time = NULL, end_time = NULL, discrete_genome = NULL,
                          keep = NULL, record_provenance = TRUE) {
  if (missing(tree_sequence) || is.null(tree_sequence)) {
    stop("First argument must be a TreeSequence instance.", call. = FALSE)
  }
  tables <- tryCatch(ts_tables(tree_sequence), error = function(e) {
    stop("First argument must be a TreeSequence instance.", call. = FALSE)
  })
  seed <- .parse_random_seed(random_seed)
  provenance_dict <- NULL
  if (record_provenance) {
    provenance_dict <- .get_provenance_dict(.parameters_dict("sim_mutations", list(
      tree_sequence = tree_sequence, rate = rate, model = model, start_time = start_time,
      end_time = end_time, discrete_genome = discrete_genome, keep = keep,
      random_seed = seed)))
  }
  rate <- rate %||% 0
  rate_map <- if (inherits(rate, "RateMap")) {
    rate
  } else if (.is_number(rate)) {
    RateMap$uniform(tables$sequence_length, as.double(rate))
  } else {
    stop("rate must be a float or a RateMap", call. = FALSE)
  }
  if (identical(tables$time_units, "uncalibrated")) {
    stop("Simulating mutations doesn't make sense when time is uncalibrated", call. = FALSE)
  }
  start_time <- if (is.null(start_time)) -.Machine$double.xmax else
    .check_number(start_time, "start_time")
  end_time <- if (is.null(end_time)) .Machine$double.xmax else .check_number(end_time, "end_time")
  if (start_time > end_time) stop("start_time must be <= end_time", call. = FALSE)
  discrete_genome <- .parse_flag(discrete_genome, default = TRUE)
  keep <- .parse_flag(keep, default = TRUE)
  model <- .mutation_model_factory(model)

  res <- .Call(msprimer_sim_mutations, tables, seed, rate_map$position, rate_map$rate,
               model$ll_spec(), discrete_genome, keep, start_time, end_time)
  model$update_state(res$model_state)
  tables <- res$tables
  if (!is.null(provenance_dict)) {
    tables <- .add_provenance(tables, .json_encode_provenance(provenance_dict))
  }
  .make_result(tables)
}

#' Simulate mutations (deprecated)
#'
#' Legacy interface simulating infinite sites binary mutations on a continuous
#' genome. Please use [sim_mutations()] in new code.
#'
#' @inheritParams sim_mutations
#' @param model `NULL` or an [InfiniteSites()] model.
#' @return A tree sequence.
#' @export
mutate <- function(tree_sequence, rate = NULL, random_seed = NULL, model = NULL,
                   keep = FALSE, start_time = NULL, end_time = NULL) {
  if (!is.null(model) && !inherits(model, "InfiniteSites")) {
    stop("This is a legacy interface which does not support general mutation models. ",
         "Please use the sim_mutations function instead.", call. = FALSE)
  }
  if (is.null(model)) model <- BinaryMutationModel()
  sim_mutations(tree_sequence, rate = rate, random_seed = random_seed, model = model,
                keep = keep, start_time = start_time, end_time = end_time,
                discrete_genome = FALSE)
}
