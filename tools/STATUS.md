# msprimeR development status (saved 2026-10-09)

All work is in a consistent, passing state; nothing is half-edited.

## Verified
- `R CMD check --as-cran`: tests, source build and all vignettes pass. Remaining
  warnings/notes are environment/repository related (missing qpdf/checkbashisms,
  site `-march` flag, no network; optional R companion packages are not on CRAN).
- testthat suite passes (parity digests in tests/testthat/fixtures).
- Full Python comparison: 66/66 scenarios bit-identical (tools/parity).
- Memory stable under stress tests (eager frees of C objects).
- Package name, native registration, options, provenance and R class tags use
  the corrected `msprimeR` spelling.
- DESCRIPTION uses the sommer maintainer identity; GPL-3 `LICENSE` is included.
  `inst/CITATION` cites canonical msprime first and identifies this R mirror;
  `citation("msprimeR")` returns exactly two entries.
- Three vignettes cover getting started, interoperability, and models/pedigrees.
- Source-built tskitR, tsdateR and tsinferR file-handoff tests passed. tsinferR's
  inference pipeline itself starts from genotype sources.

## How to resume
- Install: `R CMD INSTALL -l /tmp/msprimeR_lib . && ./cleanup`
- Regenerate docs: `Rscript -e 'roxygen2::roxygenise()' && ./cleanup`
- Re-vendor C sources: `tools/sync_c_sources.sh ../msprime ../tskit`
- Python reference (venv in /tmp may be lost on reboot): copy ../msprime to
  /tmp/msprime_ref, populate git-submodules/tskit from ../tskit (c/ and
  python/lwt_interface), set requires-python >=3.10, then
  `PYTHONPATH= CFLAGS="$(gsl-config --cflags)" LDFLAGS="$(gsl-config --libs)" pip install --no-build-isolation .`
  plus `pip install demes`; run `tools/parity/make_reference.py /tmp/parity`
  and `Rscript tools/parity/check_parity.R /tmp/parity /tmp/msprimeR_lib`.

## Remaining TODO
1. tskitR, tsdateR and tsinferR are Suggests but not all available from CRAN;
  resolve package distribution/Additional_repositories before CRAN submission.
2. Check on Windows/macOS (win-builder, R-hub), especially GSL linking.
3. Optional: port more of msprime's Python tests; the `mspms` CLI is not ported.
