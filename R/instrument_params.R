# Instrument/method-specific CentWave starting values and IPO2 search-space
# bounds, keyed by instrument x column x polarity. Edit this file to add or
# tune entries — it's kept separate from the optimization/peak-picking logic
# in R/peak_picking.R so instrument-specific numbers can be maintained
# without touching pipeline code.
#
# Falls back to IPO2's own defaults (see run_ipo_optimization() in
# R/peak_picking.R) for any (instrument, column, polarity) combo not listed
# here — fill in SynaptXS (or add other instruments) as you tune them.
#
# Each entry:
#   starting_values: ppm, peakwidth (c(min, max)), snthresh, prefilter
#     (c(k, int)), mzdiff, noise — fixed (not optimized) CentWave settings.
#   bounds: min_peakwidth, max_peakwidth, mzdiff, ppm — each c(lower, upper),
#     the IPO2 search-space range for that parameter.
#   int_threshold: absolute intensity threshold for the instrument (not used
#     by this pipeline yet — kept for the later filtering stage).

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
          mzdiff = c(-Inf, Inf), ppm = c(1, 25)
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
          mzdiff = c(-Inf, Inf), ppm = c(1, 25)
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
        min_peakwidth = c(0, 5), max_peakwidth = c(5, 45),
        mzdiff = c(-Inf, Inf), ppm = c(5, 18)
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

#' Get the `instrument` value for a group, if the sheet has that column and
#' it's actually filled in. Returns NA otherwise, so callers fall back to
#' generic behavior (e.g. IPO2's own defaults, or a relative QC threshold)
#' rather than erroring.
#'
#' @param sheet Sample sheet rows for a group.
group_instrument <- function(sheet) {
  if (!"instrument" %in% names(sheet)) {
    return(NA_character_)
  }
  values <- unique(stats::na.omit(sheet$instrument))
  if (length(values) == 0) NA_character_ else values[1]
}

#' Look up instrument/method-specific centWave config for a given
#' instrument x column x polarity combo.
#'
#' Falls through to a column-level entry (no polarity split) if a
#' polarity-specific one doesn't exist, matching how HILIC is configured
#' above for MRT. Returns NULL if nothing is configured for this combo (or
#' `instrument` is NA/unrecognized), signaling the caller to fall back to
#' generic behavior (e.g. IPO2's own defaults, or a relative QC threshold).
#'
#' @param instrument Instrument name (e.g. "MRT"), or NA.
#' @param column Chromatography column (e.g. "RP", "HILIC").
#' @param polarity Polarity ("POS" or "NEG").
#' @return A list (starting_values, bounds, int_threshold), or NULL.
get_instrument_params <- function(instrument, column, polarity) {
  if (is.na(instrument) || !instrument %in% names(INSTRUMENT_PARAMS)) {
    return(NULL)
  }

  instrument_entry <- INSTRUMENT_PARAMS[[instrument]]
  if (!column %in% names(instrument_entry)) {
    return(NULL)
  }
  column_entry <- instrument_entry[[column]]

  # Column-level entry (no polarity split), e.g. MRT/HILIC.
  if (!is.null(column_entry$starting_values)) {
    return(column_entry)
  }

  # Polarity-specific entry, e.g. MRT/RP/POS.
  if (polarity %in% names(column_entry) && !is.null(column_entry[[polarity]])) {
    return(column_entry[[polarity]])
  }

  NULL
}

#' Convert an instrument config entry (as returned by
#' `get_instrument_params()`) into the `parameter_list` shape
#' `IPO2::optimize_centwave()` expects.
#'
#' @param instrument_config One entry from `INSTRUMENT_PARAMS` (via
#'   `get_instrument_params()`).
build_ipo2_parameter_list <- function(instrument_config) {
  sv <- instrument_config$starting_values
  b <- instrument_config$bounds

  list(
    ppm = b$ppm,
    min_peakwidth = b$min_peakwidth,
    max_peakwidth = b$max_peakwidth,
    snthresh = sv$snthresh,
    prefilter_k = sv$prefilter[1],
    prefilter_int = sv$prefilter[2],
    mzdiff = b$mzdiff,
    noise = sv$noise
  )
}
