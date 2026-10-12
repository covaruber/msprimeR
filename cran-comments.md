## Test environments

* local Linux (RHEL 9), R 4.2.2, GSL 2.7

## R CMD check results

0 errors | 0 warnings (other than missing local tools: qpdf, checkbashisms) | notes:

* New submission.
* Suggests 'tskitR', which is not (yet) on CRAN; it is optional and all uses
  are guarded by `requireNamespace()`.

## Notes for CRAN

* The package bundles the msprime (GPL >= 3), tskit (MIT), kastore (MIT) and
  libavl (LGPL >= 2.1) C libraries; see inst/COPYRIGHTS. The vendored sources
  are patched (tools/sync_c_sources.sh, excluded from the build) so that
  internal assertion failures call `Rf_error()` instead of writing to stderr and
  calling `abort()`.
* SystemRequirements: GNU Scientific Library (GSL >= 2.0), located by the
  configure script via `gsl-config`, `pkg-config` or `GSL_HOME`.
