test_that("group_instrument returns NA when the column is missing or empty", {
  expect_true(is.na(group_instrument(data.frame(x = 1:3))))
  expect_true(is.na(group_instrument(data.frame(instrument = NA_character_))))
})

test_that("group_instrument returns the first non-NA instrument value", {
  sheet <- data.frame(instrument = c(NA, "MRT", "MRT"))
  expect_equal(group_instrument(sheet), "MRT")
})

test_that("get_instrument_params returns NULL for NA instrument", {
  expect_null(get_instrument_params(NA_character_, "RP", "POS"))
})

test_that("get_instrument_params matches a configured instrument/column/polarity", {
  result <- get_instrument_params("MRT", "RP", "POS")
  expect_false(is.null(result))
  expect_equal(result$starting_values$ppm, 6)
  expect_equal(result$int_threshold, 20000)
})

test_that("get_instrument_params matches case-insensitively", {
  expect_equal(get_instrument_params("mrt", "rp", "pos"), get_instrument_params("MRT", "RP", "POS"))
})

test_that("get_instrument_params falls through to a column-level entry (no polarity split)", {
  result <- get_instrument_params("MRT", "HILIC", "POS")
  expect_false(is.null(result))
  expect_equal(result$starting_values$ppm, 10)

  # Same result regardless of polarity, since HILIC has no polarity split.
  expect_equal(result, get_instrument_params("MRT", "HILIC", "NEG"))
})

test_that("get_instrument_params warns and returns NULL for an unknown instrument", {
  expect_warning(result <- get_instrument_params("Unknown", "RP", "POS"), "doesn't match any entry")
  expect_null(result)
})

test_that("get_instrument_params warns and returns NULL for an unknown column", {
  expect_warning(result <- get_instrument_params("MRT", "Unknown", "POS"), "doesn't match any configured column")
  expect_null(result)
})

test_that("get_instrument_params warns and returns NULL for an unknown polarity", {
  expect_warning(result <- get_instrument_params("MRT", "RP", "XYZ"), "doesn't match any configured polarity")
  expect_null(result)
})

test_that("get_instrument_params returns NULL silently for a deliberately-unconfigured combo", {
  expect_no_warning(result <- get_instrument_params("SynaptXS", "RP", "POS"))
  expect_null(result)
})

test_that("build_ipo2_search_space translates bounds into optimVars/lower/upper", {
  skip_if_not_installed("xcms") # build_centwave_param() calls xcms::CentWaveParam()

  instrument_config <- get_instrument_params("MRT", "RP", "POS")
  search_space <- build_ipo2_search_space(instrument_config)

  expect_setequal(search_space$optimVars, c("min_peakwidth", "max_peakwidth", "mzdiff", "ppm"))
  expect_equal(
    unname(search_space$lower[c("ppm", "mzdiff")]),
    c(1, -0.001)
  )
  expect_equal(
    unname(search_space$upper[c("ppm", "mzdiff")]),
    c(25, 0.01)
  )
})
