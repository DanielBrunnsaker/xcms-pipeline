test_that("classify_sample_type matches known prefixes case-insensitively", {
  expect_equal(classify_sample_type("sQC01"), "sQC")
  expect_equal(classify_sample_type("sqc01"), "sQC")
  expect_equal(classify_sample_type("ltQC03"), "ltQC")
  expect_equal(classify_sample_type("SolvBlank_1"), "Blank")
  expect_equal(classify_sample_type("Patient001"), "Sample")
})

test_that("extract_plate finds a Plate<N> token among hyphen-delimited parts", {
  no_plate <- extract_plate("B5W17")
  expect_equal(no_plate$value, "B5W17")
  expect_true(is.na(no_plate$plate))

  with_plate <- extract_plate("B5-Plate3")
  expect_equal(with_plate$value, "B5")
  expect_equal(with_plate$plate, "Plate3")

  only_plate <- extract_plate("Plate2")
  expect_equal(only_plate$value, "")
  expect_equal(only_plate$plate, "Plate2")
})

test_that("parse_mzml_filename parses a well-formed filename", {
  result <- parse_mzml_filename("2024-04-27_B5W17_RP_POS_sQC01_036.mzML")

  expect_equal(result$date, as.Date("2024-04-27"))
  expect_equal(result$batch, "B5W17")
  expect_equal(result$column, "RP")
  expect_equal(result$polarity, "POS")
  expect_equal(result$sample_name, "sQC01")
  expect_equal(result$sample_type, "sQC")
  expect_true(result$is_qc)
  expect_equal(result$injection_order, 36L)
  expect_true(is.na(result$plate))
  expect_equal(result$batch_plate, "B5W17")
})

test_that("parse_mzml_filename accepts an optional .raw before the extension", {
  with_raw <- parse_mzml_filename("2024-04-27_B5W17_RP_POS_sQC01_036.raw.mzML")
  without_raw <- parse_mzml_filename("2024-04-27_B5W17_RP_POS_sQC01_036.mzML")
  expect_equal(with_raw[setdiff(names(with_raw), "filename")], without_raw[setdiff(names(without_raw), "filename")])
})

test_that("parse_mzml_filename picks up a Plate token from the sample-name field", {
  result <- parse_mzml_filename("2024-04-27_B5W17_RP_POS_Plate2-sQC01_036.mzML")
  expect_equal(result$sample_name, "sQC01")
  expect_equal(result$plate, "Plate2")
  expect_equal(result$batch_plate, "B5W17_Plate2")
})

test_that("parse_mzml_filename strips a superfluous 'Plate_' marker before parsing", {
  result <- parse_mzml_filename("2024-04-27_B5W17_RP_POS_Plate_sQC01_036.mzML")
  expect_equal(result$sample_name, "sQC01")
  expect_true(is.na(result$plate))
  # Original (un-stripped) filename is preserved in the output.
  expect_equal(result$filename, "2024-04-27_B5W17_RP_POS_Plate_sQC01_036.mzML")
})

test_that("parse_mzml_filename errors on conflicting plate tokens", {
  expect_error(
    parse_mzml_filename("2024-04-27_B5-Plate2_RP_POS_Plate3-sQC01_036.mzML"),
    "Conflicting plate tokens"
  )
})

test_that("parse_mzml_filename discards ignorable extra text appended to batch", {
  result <- parse_mzml_filename("2024-06-10_B12W24_240610_plates_N29-32_RP_NEG_ltQC01_136.mzML")
  expect_equal(result$batch, "B12W24")
  expect_equal(result$column, "RP")
  expect_equal(result$polarity, "NEG")
  expect_equal(result$sample_name, "ltQC01")
  expect_equal(result$injection_order, 136L)
})

test_that("parse_mzml_filename discards ignorable extra text appended to sample name", {
  result <- parse_mzml_filename("2024-03-22_B1W12_RP_NEG_SolvBlank-P2-02_ExtraInjection_012.mzML")
  expect_equal(result$batch, "B1W12")
  expect_equal(result$sample_name, "SolvBlank-P2-02")
  expect_equal(result$sample_type, "Blank")
  expect_equal(result$injection_order, 12L)
})

test_that("parse_mzml_filename discards extra text on both batch and sample name at once", {
  result <- parse_mzml_filename(
    "2024-06-10_B12W24_240610_plates_N29-32_RP_NEG_Plate_N29-LV2006097728_203.mzML"
  )
  expect_equal(result$batch, "B12W24")
  expect_equal(result$sample_name, "N29-LV2006097728")
  expect_equal(result$injection_order, 203L)
})

test_that("parse_mzml_filename still errors on a non-matching filename (caller decides what to do)", {
  expect_error(
    parse_mzml_filename("not_a_valid_filename.mzML"),
    "does not match expected pattern"
  )
  expect_error(
    parse_mzml_filename("2024-04-27_B5W17_RP_XYZ_sQC01_036.mzML"),
    "does not match expected pattern"
  )
})

test_that("build_unparsed_row keeps the filename and salvages a leading date, blanks everything else", {
  row <- build_unparsed_row("2024-08-26_B14W34_RP_NEG_sQC_EXTRA.mzML", "some parse error")

  expect_equal(row$filename, "2024-08-26_B14W34_RP_NEG_sQC_EXTRA.mzML")
  expect_equal(row$date, as.Date("2024-08-26"))
  expect_true(row$needs_review)
  expect_equal(row$parse_error, "some parse error")
  expect_true(is.na(row$batch))
  expect_true(is.na(row$column))
  expect_true(is.na(row$polarity))
  expect_true(is.na(row$sample_name))
  expect_true(is.na(row$sample_type))
  expect_true(is.na(row$is_qc))
  expect_true(is.na(row$injection_order))
})

test_that("build_unparsed_row leaves date NA too if even that can't be found", {
  row <- build_unparsed_row("totally_unparseable.mzML", "some parse error")
  expect_true(is.na(row$date))
})

test_that("build_unparsed_row's columns match parse_mzml_filename's success-case columns", {
  ok <- parse_mzml_filename("2024-04-27_B5W17_RP_POS_sQC01_036.mzML")
  unparsed <- build_unparsed_row("bad.mzML", "err")
  expect_equal(names(ok), names(unparsed))
})

test_that("a successfully parsed row has needs_review = FALSE and parse_error = NA", {
  result <- parse_mzml_filename("2024-04-27_B5W17_RP_POS_sQC01_036.mzML")
  expect_false(result$needs_review)
  expect_true(is.na(result$parse_error))
})
