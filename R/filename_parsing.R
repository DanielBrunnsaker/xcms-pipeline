# Parsing of mzML filenames of the form:
#   DATE_BATCH_COLUMN_POLARITY_SAMPLENAME_INJORDER.raw.mzML
# e.g. 2024-03-22_B1W12_RP_POS_LV2006097259_055.raw.mzML
#      2024-03-22_B1W12_RP_POS_sQC_003.raw.mzML

MZML_FILENAME_PATTERN <- paste0(
  "^(\\d{4}-\\d{2}-\\d{2})_",  # date
  "([^_]+)_",                 # batch
  "([^_]+)_",                 # column
  "(POS|NEG)_",                # polarity
  "([^_]+)_",                  # sample name (or sQC/ltQC)
  "(\\d+)",                    # injection order
  "\\.raw\\.mzml$"
)

QC_SAMPLE_NAMES <- c("sQC", "ltQC")

#' Parse a single mzML filename into its metadata fields.
#'
#' @param filename Character scalar, just the basename (not full path).
#' @return A one-row data frame, or throws an error if the filename doesn't
#'   match the expected pattern.
parse_mzml_filename <- function(filename) {
  m <- regmatches(
    filename,
    regexec(MZML_FILENAME_PATTERN, filename, ignore.case = TRUE)
  )[[1]]

  if (length(m) == 0) {
    stop(sprintf("Filename does not match expected pattern: %s", filename), call. = FALSE)
  }

  sample_name <- m[6]
  is_qc <- sample_name %in% QC_SAMPLE_NAMES

  data.frame(
    filename = filename,
    date = as.Date(m[2]),
    batch = m[3],
    column = m[4],
    polarity = toupper(m[5]),
    sample_name = sample_name,
    is_qc = is_qc,
    qc_type = if (is_qc) sample_name else NA_character_,
    injection_order = as.integer(m[7]),
    stringsAsFactors = FALSE
  )
}

#' Recursively scan a directory for .mzML files and parse each filename.
#'
#' Fails loudly (stopping and listing every offending file) if any filename
#' does not match the expected pattern, so mis-named files are caught before
#' peak picking rather than silently mis-parsed.
#'
#' @param raw_dir Path to the project's raw data directory (may contain
#'   per-batch subfolders).
#' @return A data frame with one row per file, including a `filepath` column
#'   with the full path to each file.
scan_mzml_files <- function(raw_dir) {
  filepaths <- list.files(
    raw_dir,
    pattern = "\\.mzML$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(filepaths) == 0) {
    stop(sprintf("No .mzML files found under: %s", raw_dir), call. = FALSE)
  }

  filenames <- basename(filepaths)

  parsed <- lapply(filenames, function(f) {
    tryCatch(
      parse_mzml_filename(f),
      error = function(e) NULL
    )
  })

  failed <- filenames[vapply(parsed, is.null, logical(1))]
  if (length(failed) > 0) {
    stop(
      "The following files do not match the expected naming pattern ",
      "(DATE_BATCH_COLUMN_POLARITY_SAMPLENAME_INJORDER.raw.mzML):\n  ",
      paste(failed, collapse = "\n  "),
      call. = FALSE
    )
  }

  result <- do.call(rbind, parsed)
  result$filepath <- filepaths
  result
}
