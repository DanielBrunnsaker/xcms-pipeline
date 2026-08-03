# Parsing of mzML filenames of the form:
#   DATE_BATCH_COLUMN_POLARITY_SAMPLENAME_INJORDER.raw.mzML
# e.g. YYYY-MM-DD_BATCH_COL_POL_SAMPLENAME_INJORDER.raw.mzML
# The sample-name field may be a real sample ID, a QC type (sQC/ltQC), or a
# blank (SolvBlank). A "Plate<N>" token may also appear appended to the
# batch field or prepended to the sample-name field.

MZML_FILENAME_PATTERN <- paste0(
  "^(\\d{4}-\\d{2}-\\d{2})_",  # date
  "([^_]+)_",                 # batch
  "([^_]+)_",                 # column
  "(POS|NEG)_",                # polarity
  "([^_]+)_",                  # sample name (or sQC/ltQC)
  "(\\d+)",                    # injection order
  "\\.raw\\.mzml$"
)

# Non-sample name prefixes, matched case-insensitively against the start of
# the sample name (allowing trailing characters, e.g. a replicate number).
# Extend this list as new QC/blank types show up.
SAMPLE_TYPE_PREFIXES <- c(
  "sQC" = "sQC",
  "ltQC" = "ltQC",
  "SolvBlank" = "Blank"
)

#' Classify a sample name into a sample type via case-insensitive prefix
#' match. Anything that doesn't match a known non-sample prefix is a regular
#' "Sample".
classify_sample_type <- function(sample_name) {
  upper_name <- toupper(sample_name)
  for (prefix in names(SAMPLE_TYPE_PREFIXES)) {
    if (startsWith(upper_name, toupper(prefix))) {
      return(SAMPLE_TYPE_PREFIXES[[prefix]])
    }
  }
  "Sample"
}

PLATE_TOKEN_PATTERN <- "^Plate(\\d+)$"

#' Extract a "Plate<N>" token from a hyphen-delimited field, if present.
#'
#' The plate can show up appended to the batch field or prepended to the
#' sample-name field — this handles both by just looking for a "Plate<N>"
#' token among the hyphen-separated parts, matched case-insensitively.
#'
#' @param field Character scalar (already split out of the filename, e.g.
#'   the batch or sample-name field).
#' @return A list with `value` (the field with the plate token removed) and
#'   `plate` (normalized "Plate<N>", or NA if no plate token was found).
extract_plate <- function(field) {
  parts <- strsplit(field, "-", fixed = TRUE)[[1]]
  matches <- regmatches(parts, regexec(PLATE_TOKEN_PATTERN, parts, ignore.case = TRUE))
  digits <- vapply(matches, function(m) if (length(m) > 0) m[2] else NA_character_, character(1))

  is_plate <- !is.na(digits)
  if (!any(is_plate)) {
    return(list(value = field, plate = NA_character_))
  }

  list(
    value = paste(parts[!is_plate], collapse = "-"),
    plate = paste0("Plate", digits[is_plate][1])
  )
}

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

  batch_extract <- extract_plate(m[3])
  sample_extract <- extract_plate(m[6])

  if (!is.na(batch_extract$plate) && !is.na(sample_extract$plate) &&
      batch_extract$plate != sample_extract$plate) {
    stop(sprintf(
      "Conflicting plate tokens in filename %s: %s (in batch field) vs %s (in sample name field)",
      filename, batch_extract$plate, sample_extract$plate
    ), call. = FALSE)
  }
  plate <- if (!is.na(batch_extract$plate)) batch_extract$plate else sample_extract$plate

  batch <- batch_extract$value
  sample_name <- sample_extract$value
  sample_type <- classify_sample_type(sample_name)
  batch_plate <- if (is.na(plate)) batch else paste(batch, plate, sep = "_")

  data.frame(
    filename = filename,
    date = as.Date(m[2]),
    batch = batch,
    plate = plate,
    batch_plate = batch_plate,
    column = m[4],
    polarity = toupper(m[5]),
    sample_name = sample_name,
    sample_type = sample_type,
    is_qc = sample_type %in% c("sQC", "ltQC"),
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
#' Injection order is preferably read from each file's mzML acquisition
#' timestamp (the true run order) rather than the number encoded in the
#' filename. But if even a single file's timestamp can't be read, the whole
#' batch falls back to the filename-encoded number instead, rather than
#' mixing two different ordering sources within one sheet.
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

  n_files <- length(result$filepath)
  message(sprintf("Reading acquisition time and spectrum mode from %d file(s)...", n_files))
  pb <- utils::txtProgressBar(min = 0, max = n_files, style = 3)
  acquisition_times_raw <- numeric(n_files)
  spectrum_modes <- character(n_files)
  for (i in seq_len(n_files)) {
    acquisition_times_raw[i] <- as.numeric(get_mzml_acquisition_time(result$filepath[i]))
    spectrum_modes[i] <- detect_spectrum_mode(result$filepath[i])
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  acquisition_times <- as.POSIXct(acquisition_times_raw, origin = "1970-01-01", tz = "UTC")
  result$spectrum_mode <- spectrum_modes

  n_unknown_mode <- sum(is.na(spectrum_modes))
  if (n_unknown_mode > 0) {
    message(sprintf(
      "Spectrum mode (profile/centroid) could not be determined for %d file(s); ",
      n_unknown_mode
    ), "these will be left un-centroided during peak picking.")
  }

  if (all(!is.na(acquisition_times))) {
    message("Using mzML acquisition time for injection order.")
    result$injection_order <- rank(acquisition_times, ties.method = "first")
    result$injection_order_source <- "acquisition_time"
  } else {
    message(
      "Acquisition time unavailable for one or more files; ",
      "falling back to filename-derived injection order for all files."
    )
    result$injection_order_source <- "filename"
  }

  result
}
