# Standalone helper: generate a sample sheet from raw mzML filenames.
#
# This is not a pipeline stage — it's an optional convenience for deriving a
# starting sample sheet when you don't already have one. The sheet always
# needs manual curation afterwards regardless of how it was produced.
#
# Usage:
#   Rscript scripts/generate_sample_sheet.R <folder>
#
# <folder> is any folder containing your .mzML files (flat or nested, both
# fine — batch is parsed from the filename, not from folder structure). It
# does not need to live inside this repo. Writes the sheet to
# <folder>/metadata/sample_sheet.xlsx. Review and edit that file (fill in
# `sample_group`, add `notes`) before running run_peak_picking.R.

source("R/acquisition_time.R")
source("R/spectrum_mode.R")
source("R/filename_parsing.R")
source("R/sample_sheet.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scripts/generate_sample_sheet.R <folder>", call. = FALSE)
}
folder <- args[1]

if (!dir.exists(folder)) {
  stop(sprintf("Folder does not exist: %s", folder), call. = FALSE)
}

out_path <- file.path(folder, "metadata", "sample_sheet.xlsx")

generate_sample_sheet(folder, out_path)

message(
  "\nNext step: open the sample sheet, review the parsed metadata, and fill ",
  "in `sample_group` (and any `notes`) before running run_peak_picking.R."
)
