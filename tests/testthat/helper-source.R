# Sources the R/ files needed for tests, using the project root captured
# by tests/testthat.R (options("xcms_pipeline.project_root")) rather than
# a relative path -- testthat::test_dir() changes the working directory
# internally, so a plain "R/..." path wouldn't resolve here.
#
# R/parallel.R, R/peak_picking.R, and R/spectrum_mode.R each have
# Bioconductor-requiring session setup (library(xcms)/library(MSnbase), the
# IPO2 monkeypatch, BiocParallel::register()) deliberately positioned LAST in
# the file, after all function definitions -- so source()ing them here still
# defines their pure functions even without xcms/MSnbase/IPO2 installed, or
# (for parallel.R specifically) in a sandbox that can't bind a SnowParam
# socket; only that trailing setup code errors, caught below.

project_root <- getOption("xcms_pipeline.project_root", ".")

r_path <- function(file) file.path(project_root, "R", file)

source_tolerant <- function(path) {
  tryCatch(
    source(path),
    error = function(e) {
      message(sprintf("(%s not fully sourced: %s -- functions defined above the error are still available)",
                       path, conditionMessage(e)))
    }
  )
}

source(r_path("acquisition_time.R"))
source(r_path("filename_parsing.R"))
source(r_path("sample_sheet.R"))
source(r_path("instrument_params.R"))
source_tolerant(r_path("parallel.R"))
source(r_path("qc_quality.R"))
source_tolerant(r_path("spectrum_mode.R"))
source_tolerant(r_path("peak_picking.R"))
source_tolerant(r_path("retgroup_optimization.R"))
