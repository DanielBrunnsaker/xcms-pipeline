# Parsing of mzML filenames:
#   DATE_BATCH_COLUMN_POLARITY_SAMPLENAME_INJORDER[.raw].mzML
# Notes:
# - ".raw" before the extension is optional (present on some exports, not others)
# - sample-name may be a real sample ID, a QC type (sQC/ltQC), or a blank (SolvBlank)
# - "Plate<N>" may appear appended to batch, or prepended to sample-name
# - a literal "Plate_" marker is stripped before parsing (see SUPERFLUOUS_FILENAME_MARKERS)

MZML_FILENAME_PATTERN <- paste0(
  "^(\\d{4}-\\d{2}-\\d{2})_",  # date
  "([^_]+)_",                 # batch
  "([^_]+)_",                 # column
  "(POS|NEG)_",                # polarity
  "([^_]+)_",                  # sample name (or sQC/ltQC)
  "(\\d+)",                    # injection order
  "(?:\\.raw)?\\.mzml$"        # ".raw" before the extension is optional
)

# TODO: revisit — assumed to carry no information (e.g. whether the ID
# after it, like "N10", ever matters), not a confirmed-safe assumption.
SUPERFLUOUS_FILENAME_MARKERS <- c("Plate_")

# Non-sample-name prefixes, matched case-insensitively against the start of
# the sample name (allows trailing chars, e.g. a replicate number). Extend
# as new QC/blank types show up.
SAMPLE_TYPE_PREFIXES <- c(
  "sQC" = "sQC",
  "ltQC" = "ltQC",
  "SolvBlank" = "Blank"
)

#' Classify a sample name via case-insensitive prefix match against
#' SAMPLE_TYPE_PREFIXES; anything unmatched is a regular "Sample".
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

#' Extract a "Plate<N>" token from a hyphen-delimited field (batch or
#' sample-name), if present.
#'
#' @param field Character scalar already split out of the filename.
#' @return list(value = field with token removed, plate = "Plate<N>" or NA).
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
  parse_target <- filename
  for (marker in SUPERFLUOUS_FILENAME_MARKERS) {
    parse_target <- sub(marker, "", parse_target, fixed = TRUE)
  }

  m <- regmatches(
    parse_target,
    regexec(MZML_FILENAME_PATTERN, parse_target, ignore.case = TRUE)
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
#' - Fails loudly, listing every offending file, rather than silently
#'   mis-parsing a mis-named one.
#' - injection_order prefers each file's mzML acquisition timestamp (true
#'   run order) over the filename-encoded number; if even one file's
#'   timestamp is unreadable, the whole batch falls back to the
#'   filename-encoded number instead of mixing two ordering sources.
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
