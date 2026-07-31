# Stage 1: generate the sample sheet for a project from raw mzML filenames.
#
# Usage:
#   Rscript scripts/01_generate_sample_sheet.R <project_name>
#
# Reads data/<project>/raw/**/*.mzML, derives metadata from filenames, and
# writes data/<project>/metadata/sample_sheet.xlsx. Review and edit that file
# (fill in `sample_group`, add `notes`) before running 02_run_peak_picking.R.

source("R/filename_parsing.R")
source("R/sample_sheet.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/01_generate_sample_sheet.R <project_name>", call. = FALSE)
}
project <- args[1]

raw_dir <- file.path("data", project, "raw")
out_path <- file.path("data", project, "metadata", "sample_sheet.xlsx")

generate_sample_sheet(raw_dir, out_path)

message(
  "\nNext step: open the sample sheet, review the parsed metadata, and fill ",
  "in `sample_group` (and any `notes`) before running 02_run_peak_picking.R."
)
