# Test runner for the pure/deterministic logic (filename parsing, sample
# sheet handling, instrument param lookup, IPO subset selection). Doesn't
# need Docker or Bioconductor -- just testthat:
#   install.packages("testthat")
#
# Run from the repo root:
#   Rscript tests/testthat.R
#
# Anything that actually calls xcms/MSnbase/IPO2/IPO (real peak picking,
# alignment, optimization) is out of scope here -- covered by an actual
# pipeline run in Docker instead, not a fast local unit-test loop.

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop('testthat not installed. Run: install.packages("testthat")', call. = FALSE)
}

# test_dir() changes the working directory internally, so the repo root is
# captured here (while it's still the caller's CWD) for helper-source.R to
# use instead of relying on relative paths at source() time.
options(xcms_pipeline.project_root = normalizePath("."))

testthat::test_dir("tests/testthat", reporter = "summary")
