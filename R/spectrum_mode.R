# Detecting whether an mzML file is profile or centroid mode, and
# centroiding profile data before peak picking — xcms's centWave works much
# better on centroided data than profile data.
#
# MSnbase::pickPeaks() on an OnDiskMSnExp doesn't run immediately — it
# registers a lazy ProcessingStep (the function stored just as the name
# "pickPeaks") to be applied later when spectra are actually read. That
# name only resolves if MSnbase is attached, not just namespace-referenced
# via MSnbase:: — so, unlike the rest of this codebase, MSnbase needs to
# actually be library()'d for centroiding to work.
library(MSnbase)

#' Detect whether an mzML file is profile or centroid mode.
#'
#' mzML files declare this per-spectrum via a "profile spectrum" or
#' "centroid spectrum" cvParam, appearing further into the file than the
#' header alone — so this reads a larger chunk from the start of the file
#' (enough to reach the first spectrum's cvParams) and just checks for a
#' literal mention of "profile" or "centroid" (case-insensitive).
#'
#' @param filepath Path to an .mzML file.
#' @param n_bytes How many bytes to read from the start of the file.
#' @return "profile", "centroid", or NA if neither (or both) are found.
detect_spectrum_mode <- function(filepath, n_bytes = 200000) {
  tryCatch({
    con <- file(filepath, open = "rb")
    on.exit(close(con))
    raw <- readBin(con, what = "raw", n = n_bytes)
    text <- tolower(rawToChar(raw))

    has_profile <- grepl("profile", text, fixed = TRUE)
    has_centroid <- grepl("centroid", text, fixed = TRUE)

    if (has_centroid && !has_profile) {
      "centroid"
    } else if (has_profile && !has_centroid) {
      "profile"
    } else {
      NA_character_
    }
  }, error = function(e) NA_character_)
}

#' Get a sheet's `spectrum_mode` values, defaulting to NA for every row if
#' the column doesn't exist (e.g. a hand-made sheet, or one generated before
#' this detection existed) — callers treat NA as "unknown, don't centroid"
#' rather than guessing.
#'
#' @param sheet Sample sheet rows (or a subset).
get_spectrum_modes <- function(sheet) {
  if ("spectrum_mode" %in% names(sheet)) sheet$spectrum_mode else rep(NA_character_, nrow(sheet))
}

#' Read raw MS data for a set of files, centroiding any detected as profile
#' mode via `MSnbase::pickPeaks()` before returning.
#'
#' Warns (rather than silently guessing) if a read mixes profile and
#' centroid files, since that's unusual and applying centroiding
#' indiscriminately to already-centroided spectra can distort them.
#'
#' @param filepaths Character vector of mzML file paths.
#' @param spectrum_modes Character vector ("profile"/"centroid"/NA), one per
#'   file, in the same order as `filepaths`.
#' @return An OnDiskMSnExp, centroided if needed.
read_raw_data <- function(filepaths, spectrum_modes) {
  raw_data <- MSnbase::readMSData(filepaths, mode = "onDisk")

  modes_present <- unique(stats::na.omit(spectrum_modes))
  if (length(modes_present) > 1) {
    warning(
      "Mixed profile/centroid files in the same read (",
      paste(modes_present, collapse = ", "),
      "); centroiding will be applied to all of them, which isn't correct ",
      "for files that are already centroided.",
      call. = FALSE
    )
  }

  if ("profile" %in% modes_present) {
    message("Profile-mode data detected; centroiding with MSnbase::pickPeaks()...")
    raw_data <- MSnbase::pickPeaks(raw_data)
  }

  raw_data
}
