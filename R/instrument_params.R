# Instrument/method-specific CentWave starting values and IPO2 search-space
# bounds, keyed by instrument x column x polarity.
# - Edit this file to add/tune entries; kept separate from R/peak_picking.R
#   so instrument-specific numbers don't mix with pipeline logic.
# - Falls back to default_ipo2_search_space() (below) for any combo not
#   listed here — fill in SynaptXS (or other instruments) as you tune them.
#
# Each entry:
# - starting_values: ppm, peakwidth (c(min, max)), snthresh,
#   prefilter (c(k, int)), mzdiff, noise — fixed (not optimized) settings.
# - bounds: min_peakwidth, max_peakwidth, mzdiff, ppm — each
#   c(lower, upper), the IPO2 search-space range for that parameter.
# - int_threshold: absolute intensity floor for the instrument, used as
#   check_qc_quality()'s TIC threshold when available (R/qc_quality.R).

INSTRUMENT_PARAMS <- list(
  MRT = list(
    RP = list(
      POS = list(
        starting_values = list(
          ppm = 6, peakwidth = c(1.7, 22), snthresh = 10,
          prefilter = c(5, 45000), mzdiff = 0.0017, noise = 20000
        ),
        bounds = list(
          min_peakwidth = c(0.5, 5), max_peakwidth = c(10, 45),
          mzdiff = c(-0.001, 0.01), ppm = c(1, 25)
        ),
        int_threshold = 20000
      ),
      NEG = list(
        starting_values = list(
          ppm = 6, peakwidth = c(1.7, 22), snthresh = 10,
          prefilter = c(5, 45000), mzdiff = 0.0017, noise = 20000
        ),
        bounds = list(
          min_peakwidth = c(0.5, 5), max_peakwidth = c(10, 45),
          mzdiff = c(-0.001, 0.01), ppm = c(1, 25)
        ),
        int_threshold = 20000
      )
    ),
    HILIC = list(
      starting_values = list(
        ppm = 10, peakwidth = c(2.5, 30), snthresh = 10,
        prefilter = c(5, 15000), mzdiff = -0.001, noise = 10000
      ),
      bounds = list(
        min_peakwidth = c(0.5, 4), max_peakwidth = c(5, 45),
        # mzdiff lower bound widened slightly below the -0.001 starting
        # value -- IPO2::optimXCMS() requires the start to sit strictly
        # inside its bounds (start > lower), not just within them.
        mzdiff = c(-0.0015, 0.01), ppm = c(5, 18)
      ),
      int_threshold = 10000
    )
  ),
  SynaptXS = list(
    RP = list(
      POS = NULL,
      NEG = NULL
    ),
    HILIC = list(
      POS = NULL,
      NEG = NULL
    )
  )
)

#' Get the `instrument` value for a group (NA if the column is missing or
#' unfilled, so callers fall back to generic behavior instead of erroring).
#'
#' @param sheet Sample sheet rows for a group.
group_instrument <- function(sheet) {
  if (!"instrument" %in% names(sheet)) {
    return(NA_character_)
  }
  values <- unique(stats::na.omit(sheet$instrument))
  if (length(values) == 0) NA_character_ else values[1]
}

#' Look up instrument/method-specific centWave config for an
#' instrument x column x polarity combo. Matching is case-insensitive
#' (e.g. "mrt" matches "MRT").
#' - Falls through to a column-level entry (no polarity split) if no
#'   polarity-specific one exists (matches how HILIC is set up for MRT).
#' - Returns NULL if nothing is configured (or `instrument` is NA) —
#'   caller falls back to generic behavior.
#' - Warns (rather than silently falling back) if `instrument` or `column`
#'   was actually supplied but doesn't match anything known — likely a
#'   typo, not an intentional gap. Does *not* warn when the structure
#'   matches but the specific entry is deliberately `NULL` (e.g. the
#'   `SynaptXS` placeholders) — that's an expected, not-yet-tuned combo.
#'
#' @param instrument Instrument name (e.g. "MRT"), or NA.
#' @param column Chromatography column (e.g. "RP", "HILIC").
#' @param polarity Polarity ("POS" or "NEG").
#' @return A list (starting_values, bounds, int_threshold), or NULL.
get_instrument_params <- function(instrument, column, polarity) {
  if (is.na(instrument)) {
    return(NULL)
  }

  instrument_match <- names(INSTRUMENT_PARAMS)[toupper(names(INSTRUMENT_PARAMS)) == toupper(instrument)]
  if (length(instrument_match) == 0) {
    warning(sprintf(
      "instrument '%s' doesn't match any entry in INSTRUMENT_PARAMS (%s); using IPO2 defaults.",
      instrument, paste(names(INSTRUMENT_PARAMS), collapse = ", ")
    ), call. = FALSE)
    return(NULL)
  }
  instrument_entry <- INSTRUMENT_PARAMS[[instrument_match[1]]]

  column_match <- names(instrument_entry)[toupper(names(instrument_entry)) == toupper(column)]
  if (length(column_match) == 0) {
    warning(sprintf(
      "column '%s' doesn't match any configured column for instrument '%s' (%s); using IPO2 defaults.",
      column, instrument, paste(names(instrument_entry), collapse = ", ")
    ), call. = FALSE)
    return(NULL)
  }
  column_entry <- instrument_entry[[column_match[1]]]

  # Column-level entry (no polarity split), e.g. MRT/HILIC.
  if (!is.null(column_entry$starting_values)) {
    return(column_entry)
  }

  # Polarity-specific entry, e.g. MRT/RP/POS. May legitimately be NULL
  # (not yet configured) — that's not a mismatch, so no warning here.
  polarity_match <- names(column_entry)[toupper(names(column_entry)) == toupper(polarity)]
  if (length(polarity_match) > 0) {
    return(column_entry[[polarity_match[1]]])
  }

  warning(sprintf(
    "polarity '%s' doesn't match any configured polarity for %s/%s (%s); using IPO2 defaults.",
    polarity, instrument, column, paste(names(column_entry), collapse = ", ")
  ), call. = FALSE)
  NULL
}

#' Flatten a CentWaveParam into a single-row data.frame -- the inverse of
#' `build_centwave_param()`, used for human-readable summaries of which
#' parameters were actually used where (e.g. one row per batch, see
#' scripts/run_peak_picking.R's "batch" scope).
#'
#' @param centwave_param An xcms::CentWaveParam.
#' @return A one-row data.frame.
centwave_param_to_row <- function(centwave_param) {
  data.frame(
    ppm = centwave_param@ppm,
    min_peakwidth = centwave_param@peakwidth[1],
    max_peakwidth = centwave_param@peakwidth[2],
    snthresh = centwave_param@snthresh,
    prefilter_k = centwave_param@prefilter[1],
    prefilter_int = centwave_param@prefilter[2],
    mzdiff = centwave_param@mzdiff,
    noise = centwave_param@noise,
    stringsAsFactors = FALSE
  )
}

#' Build a starting-point CentWaveParam from an instrument config's
#' starting_values. Fields not covered by starting_values (mzCenterFun,
#' integrate, fitgauss, verboseColumns, roiList, firstBaselineCheck,
#' roiScales) use xcms::CentWaveParam()'s own defaults.
#'
#' @param starting_values `starting_values` list from an INSTRUMENT_PARAMS
#'   entry (or a generic fallback for unconfigured instruments).
build_centwave_param <- function(starting_values) {
  sv <- starting_values
  xcms::CentWaveParam(
    ppm = sv$ppm,
    peakwidth = sv$peakwidth,
    snthresh = sv$snthresh,
    prefilter = sv$prefilter,
    mzdiff = sv$mzdiff,
    noise = sv$noise
  )
}

#' Convert an instrument config entry (as returned by
#' `get_instrument_params()`) into the shape `IPO2::optimXCMS()` expects:
#' a starting-point CentWaveParam plus which fields to search over and
#' their bounds.
#'
#' Bounds field names (`min_peakwidth`, `max_peakwidth`, `mzdiff`, `ppm`)
#' already match `optimXCMS()`'s `optimVars` naming directly, so this is a
#' near 1:1 translation of R/instrument_params.R's existing bounds shape.
#'
#' @param instrument_config One entry from `INSTRUMENT_PARAMS` (via
#'   `get_instrument_params()`).
#' @return A list(cwParam, optimVars, lower, upper).
build_ipo2_search_space <- function(instrument_config) {
  b <- instrument_config$bounds

  list(
    cwParam = build_centwave_param(instrument_config$starting_values),
    optimVars = names(b),
    lower = vapply(b, `[`, numeric(1), 1),
    upper = vapply(b, `[`, numeric(1), 2)
  )
}

#' Generic fallback search space for column x polarity groups with no
#' instrument-specific config — IPO2 (Carl Brunius's) has no built-in
#' "suggest defaults" the way the earlier IPO2 (wmoldham's) did, so this
#' mirrors that package's suggested starting ranges by hand.
default_ipo2_search_space <- function() {
  build_ipo2_search_space(list(
    starting_values = list(
      ppm = 25, peakwidth = c(20, 50), snthresh = 100,
      prefilter = c(3, 10000), mzdiff = -0.001, noise = 0
    ),
    bounds = list(
      min_peakwidth = c(12, 28), max_peakwidth = c(35, 65),
      # Lower bound widened slightly below the -0.001 starting value -- see
      # the same note on the MRT/HILIC entry above.
      mzdiff = c(-0.0015, 0.01), ppm = c(17, 32)
    )
  ))
}
