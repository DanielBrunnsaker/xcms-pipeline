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

test_that("parse_mzml_filename errors on a non-matching filename", {
  expect_error(
    parse_mzml_filename("not_a_valid_filename.mzML"),
    "does not match expected pattern"
  )
  expect_error(
    parse_mzml_filename("2024-04-27_B5W17_RP_XYZ_sQC01_036.mzML"),
    "does not match expected pattern"
  )
})
