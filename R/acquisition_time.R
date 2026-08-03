# Reading the acquisition start time directly from mzML XML, without parsing
# the full (potentially large) document.

#' Read the acquisition start time from an mzML file's <run> element.
#'
#' The `startTimeStamp` attribute sits near the top of the file, well before
#' the (potentially huge) spectrum data, so this only reads a small chunk of
#' the file rather than the whole thing.
#'
#' @param filepath Path to an .mzML file.
#' @param n_bytes How many bytes to read from the start of the file when
#'   looking for the <run> element.
#' @return A POSIXct, or NA if the timestamp couldn't be found or parsed.
get_mzml_acquisition_time <- function(filepath, n_bytes = 20000) {
  tryCatch({
    con <- file(filepath, open = "rb")
    on.exit(close(con))
    raw <- readBin(con, what = "raw", n = n_bytes)
    text <- rawToChar(raw)

    m <- regmatches(text, regexpr('startTimeStamp="[^"]+"', text))
    if (length(m) == 0) {
      return(NA)
    }

    timestamp_str <- sub('startTimeStamp="([^"]+)"', "\\1", m)
    as.POSIXct(timestamp_str, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
  }, error = function(e) NA)
}
