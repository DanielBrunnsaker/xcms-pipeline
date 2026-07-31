# Assembly of the sample sheet written for manual review.

#' Build the sample sheet data frame from parsed mzML filename metadata.
#'
#' Adds columns that can't be derived from the filename (left blank for the
#' user to fill in during manual review) and sorts rows for easy checking.
#'
#' @param parsed Data frame as returned by `scan_mzml_files()`.
#' @return A data frame ready to be written to Excel.
build_sample_sheet <- function(parsed) {
  parsed$sample_group <- NA_character_
  parsed$notes <- NA_character_

  parsed <- parsed[order(parsed$column, parsed$polarity, parsed$batch, parsed$injection_order), ]

  parsed[, c(
    "filepath", "filename", "date", "batch", "column", "polarity",
    "sample_name", "is_qc", "qc_type", "injection_order",
    "sample_group", "notes"
  )]
}

#' Generate and write the sample sheet for a project.
#'
#' @param raw_dir Path to the project's raw data directory.
#' @param out_path Path to write the sample sheet Excel file to.
generate_sample_sheet <- function(raw_dir, out_path) {
  parsed <- scan_mzml_files(raw_dir)
  sheet <- build_sample_sheet(parsed)

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(sheet, out_path)

  message(sprintf("Wrote sample sheet with %d rows to: %s", nrow(sheet), out_path))
  invisible(sheet)
}
