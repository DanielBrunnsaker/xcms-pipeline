#' Detect whether an mzML file is profile or centroid mode.
#'
#' mzML declares this per-spectrum via a "profile spectrum"/"centroid
#' spectrum" cvParam, past the header — reads a larger chunk from the file
#' start (enough to reach the first spectrum's cvParams) and checks for a
#' literal "profile"/"centroid" mention (case-insensitive).
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
#' Also strips completely empty spectra (zero data points -- not "few
#' peaks", nothing at all) via `MSnbase::filterEmptySpectra()`. A trailing
#' empty spectrum crashes `adjustRtime(ObiwarpParam())`'s profile-matrix
#' construction with "subscript out of bounds" -- a real, still-unfixed
#' xcms bug (github.com/sneumann/xcms/issues/366) in `.insertColumn()`'s
#' handling of that case, confirmed against xcms's current source, not a
#' data-quality problem with the file otherwise. An xcms maintainer's
#' recommended workaround for exactly this crash, confirmed by another user
#' in the same thread, is filtering empty spectra out before peak
#' picking/alignment -- applied here, once, for every raw-data read in the
#' pipeline (both peak-picking and QC checking share this function), rather
#' than requiring a caller-by-caller fix.
#'
#' @param filepaths Character vector of mzML file paths.
#' @param spectrum_modes Character vector ("profile"/"centroid"/NA), one per
#'   file, in the same order as `filepaths`.
#' @return An OnDiskMSnExp, centroided if needed, with empty spectra removed.
read_raw_data <- function(filepaths, spectrum_modes) {
  raw_data <- MSnbase::readMSData(filepaths, mode = "onDisk")

  n_before <- length(raw_data)
  raw_data <- MSnbase::filterEmptySpectra(raw_data)
  n_removed <- n_before - length(raw_data)
  if (n_removed > 0) {
    message(sprintf(
      "Removed %d completely empty spectrum/spectra (of %d) -- known adjustRtime()/obiwarp crash trigger, see read_raw_data().",
      n_removed, n_before
    ))
  }

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

# MSnbase::pickPeaks() on an OnDiskMSnExp doesn't run immediately -- it
# registers a lazy ProcessingStep (function stored as the name "pickPeaks")
# to apply later when spectra are read. That name only resolves if MSnbase
# is attached, not just namespace-referenced. Deliberately last in this
# file (not first): has no effect on the functions defined above, and
# keeps them sourceable/testable without MSnbase installed.
library(MSnbase)
