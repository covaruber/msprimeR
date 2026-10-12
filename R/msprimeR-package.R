#' msprimeR: coalescent simulations with msprime, natively in R
#'
#' An independent native R mirror of the 'msprime' package. All simulation work
#' is performed by the bundled original 'msprime' C library; this package
#' provides an R interface for demography, ancestry and mutation models, rate
#' maps, pedigrees and likelihoods. This is not an official upstream msprime
#' release. Cite the canonical msprime 1.0 publication: Baumdicker et al. (2022),
#' *Efficient ancestry and mutation simulation with msprime 1.0*, Genetics
#' 220(3), iyab229, doi:10.1093/genetics/iyab229. See `citation("msprimeR")`.
#'
#' Population, individual and node IDs are **zero-based**, exactly as in the
#' tree sequence tables (and in Python's msprime and the 'tskitR' package).
#'
#' @keywords internal
#' @useDynLib msprimeR, .registration = TRUE
#' @import R6
#' @importFrom stats runif setNames
#' @importFrom utils head tail
"_PACKAGE"
