# Assembly of the sample sheet written for manual review.

#' Copy a file to a timestamped backup path before it gets overwritten, if
#' it exists. Guards against silently clobbering manual edits (sample_group,
#' notes, QC review, ...) when a script re-writes the same sheet path.
#'
#' @param path Path to the file that's about to be overwritten.
backup_file <- function(path) {
  if (!file.exists(path)) {
    return(invisible(NULL))
  }

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  backup_path <- sub("(\\.[^.]+)$", paste0("_backup_", timestamp, "\\1"), path)
  file.copy(path, backup_path)
  message(sprintf("Backed up existing file to: %s", backup_path))
  invisible(backup_path)
}

#' Disambiguate sample names repeated within a batch (e.g. two "sQC01"
#' injections not recorded with a plate in the filename).
#' - Infers a synthetic plate from acquisition order: 1st occurrence ->
#'   "Plate1-sQC01", 2nd -> "Plate2-sQC01", etc. Uses the real `plate`
#'   instead if already parsed. Unique names are left as-is.
#' - A heuristic, not ground truth -- kept in a separate `sample_label`
#'   column (not overwriting `sample_name`/`plate`) so a wrong guess is easy
#'   to spot and fix during review.
#'
#' @param parsed Data frame as returned by `scan_mzml_files()`.
#' @return The same data frame with a new `sample_label` column.
disambiguate_sample_names <- function(parsed) {
  parsed$sample_label <- parsed$sample_name

  for (idx in split(seq_len(nrow(parsed)), parsed$batch)) {
    batch_names <- parsed$sample_name[idx]
    dup_names <- unique(batch_names[duplicated(batch_names)])

    for (name in dup_names) {
      name_idx <- idx[batch_names == name]
      name_idx <- name_idx[order(parsed$injection_order[name_idx])]

      for (i in seq_along(name_idx)) {
        row_i <- name_idx[i]
        plate_label <- if (!is.na(parsed$plate[row_i])) parsed$plate[row_i] else paste0("Plate", i)
        parsed$sample_label[row_i] <- paste0(plate_label, "-", parsed$sample_name[row_i])
      }
    }
  }

  parsed
}

#' Build the sample sheet data frame from parsed mzML filename metadata.
#'
#' Adds columns that can't be derived from the filename (left blank for the
#' user to fill in during manual review) and sorts rows for easy checking.
#'
#' `include` defaults to FALSE (not the usual TRUE) for `needs_review` rows
#' -- a file whose name couldn't be parsed at all has no reliable
#' batch/column/polarity/sample_type, so it shouldn't silently enter QC
#' checking or peak picking until someone's filled those in by hand and
#' flipped `include` back on.
#'
#' @param parsed Data frame as returned by `scan_mzml_files()`.
#' @return A data frame ready to be written to Excel.
build_sample_sheet <- function(parsed) {
  parsed <- disambiguate_sample_names(parsed)
  parsed$instrument <- NA_character_
  parsed$sample_group <- NA_character_
  parsed$notes <- NA_character_
  parsed$include <- !parsed$needs_review

  parsed <- parsed[order(parsed$column, parsed$polarity, parsed$batch_plate, parsed$injection_order), ]

  parsed[, c(
    "filepath", "filename", "needs_review", "parse_error", "date", "batch",
    "plate", "batch_plate", "column", "polarity", "sample_name",
    "sample_label", "sample_type", "is_qc", "injection_order",
    "injection_order_source", "spectrum_mode", "instrument", "sample_group",
    "notes", "include"
  )]
}

#' Get a sheet's `include` values, defaulting to TRUE for every row if the
#' column doesn't exist (e.g. a hand-made sheet, or one generated before
#' this existed) or is blank for a given row — inclusion is opt-out, not
#' opt-in, so absence never silently excludes anything.
#'
#' @param sheet Sample sheet rows (or a subset).
get_included <- function(sheet) {
  if (!"include" %in% names(sheet)) {
    return(rep(TRUE, nrow(sheet)))
  }
  included <- as.logical(sheet$include)
  included[is.na(included)] <- TRUE
  included
}

#' Generate and write the sample sheet for a project.
#'
#' @param raw_dir Path to the project's raw data directory.
#' @param out_path Path to write the sample sheet Excel file to.
generate_sample_sheet <- function(raw_dir, out_path) {
  parsed <- scan_mzml_files(raw_dir)
  sheet <- build_sample_sheet(parsed)

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  backup_file(out_path)
  writexl::write_xlsx(sheet, out_path)

  message(sprintf("Wrote sample sheet with %d rows to: %s", nrow(sheet), out_path))
  invisible(sheet)
}

#' Validate that a sample sheet has what the peak-picking pipeline needs.
#'
#' The sheet may have been produced by `generate_sample_sheet()` or supplied
#' by hand, so this checks the columns the pipeline actually depends on
#' rather than assuming any particular origin.
#'
#' @param sheet Data frame read from a sample sheet Excel file.
validate_sample_sheet <- function(sheet) {
  required_cols <- c(
    "filepath", "column", "polarity", "sample_name", "sample_label",
    "sample_type", "is_qc", "injection_order"
  )
  missing_cols <- setdiff(required_cols, names(sheet))
  if (length(missing_cols) > 0) {
    stop(
      "Sample sheet is missing required column(s): ", paste(missing_cols, collapse = ", "),
      ". Required columns: ", paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.character(sheet$filepath)) {
    stop("Sample sheet column `filepath` must be character.", call. = FALSE)
  }
  if (!is.character(sheet$sample_name)) {
    stop("Sample sheet column `sample_name` must be character.", call. = FALSE)
  }
  if (!is.character(sheet$sample_label)) {
    stop("Sample sheet column `sample_label` must be character.", call. = FALSE)
  }
  if (!is.character(sheet$sample_type)) {
    stop("Sample sheet column `sample_type` must be character.", call. = FALSE)
  }
  if (!is.logical(sheet$is_qc)) {
    stop("Sample sheet column `is_qc` must be logical (TRUE/FALSE).", call. = FALSE)
  }
  if (!is.numeric(sheet$injection_order)) {
    stop("Sample sheet column `injection_order` must be numeric.", call. = FALSE)
  }

  missing_files <- sheet$filepath[!file.exists(sheet$filepath)]
  if (length(missing_files) > 0) {
    stop(
      "Sample sheet references file(s) that don't exist on disk:\n  ",
      paste(missing_files, collapse = "\n  "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
