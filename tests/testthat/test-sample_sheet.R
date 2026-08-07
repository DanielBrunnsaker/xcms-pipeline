test_that("get_included defaults to TRUE when the column is missing", {
  sheet <- data.frame(x = 1:3)
  expect_equal(get_included(sheet), c(TRUE, TRUE, TRUE))
})

test_that("get_included treats NA as TRUE (opt-out, not opt-in)", {
  sheet <- data.frame(include = c(TRUE, FALSE, NA))
  expect_equal(get_included(sheet), c(TRUE, FALSE, TRUE))
})

test_that("disambiguate_sample_names leaves unique names alone", {
  parsed <- data.frame(
    batch = c("B1", "B1"),
    sample_name = c("Sample1", "Sample2"),
    plate = NA_character_,
    injection_order = c(1, 2),
    stringsAsFactors = FALSE
  )
  result <- disambiguate_sample_names(parsed)
  expect_equal(result$sample_label, c("Sample1", "Sample2"))
})

test_that("disambiguate_sample_names labels duplicates by chronological plate order", {
  parsed <- data.frame(
    batch = c("B1", "B1", "B1"),
    sample_name = c("sQC01", "sQC01", "sQC01"),
    plate = NA_character_,
    injection_order = c(30, 10, 20),
    stringsAsFactors = FALSE
  )
  result <- disambiguate_sample_names(parsed)
  # Row 2 (injection_order 10) is chronologically first -> Plate1, row 3
  # (order 20) second -> Plate2, row 1 (order 30) third -> Plate3.
  expect_equal(result$sample_label[2], "Plate1-sQC01")
  expect_equal(result$sample_label[3], "Plate2-sQC01")
  expect_equal(result$sample_label[1], "Plate3-sQC01")
})

test_that("disambiguate_sample_names uses a real plate over a synthetic one", {
  parsed <- data.frame(
    batch = c("B1", "B1"),
    sample_name = c("sQC01", "sQC01"),
    plate = c("Plate9", NA),
    injection_order = c(1, 2),
    stringsAsFactors = FALSE
  )
  result <- disambiguate_sample_names(parsed)
  expect_equal(result$sample_label[1], "Plate9-sQC01")
})

test_that("build_sample_sheet defaults include to FALSE for needs_review rows, TRUE otherwise", {
  parsed <- data.frame(
    filepath = c("f1.mzML", "f2.mzML"), filename = c("f1.mzML", "f2.mzML"),
    date = as.Date(c("2024-01-01", "2024-01-02")),
    batch = c("B1", NA), plate = NA_character_, batch_plate = c("B1", NA),
    column = c("RP", NA), polarity = c("POS", NA),
    sample_name = c("s1", NA), sample_type = c("Sample", NA),
    is_qc = c(FALSE, NA), injection_order = c(1L, NA),
    injection_order_source = "filename", spectrum_mode = NA_character_,
    needs_review = c(FALSE, TRUE), parse_error = c(NA_character_, "bad filename"),
    stringsAsFactors = FALSE
  )
  sheet <- build_sample_sheet(parsed)
  expect_equal(sheet$include[sheet$filename == "f1.mzML"], TRUE)
  expect_equal(sheet$include[sheet$filename == "f2.mzML"], FALSE)
})

test_that("validate_sample_sheet errors on a missing required column", {
  sheet <- data.frame(filepath = character(0))
  expect_error(validate_sample_sheet(sheet), "missing required column")
})

test_that("validate_sample_sheet errors on wrong column types", {
  f <- tempfile(fileext = ".mzML")
  file.create(f)
  on.exit(unlink(f))

  base <- data.frame(
    filepath = f, column = "RP", polarity = "POS", sample_name = "s1",
    sample_label = "s1", sample_type = "Sample", is_qc = FALSE,
    injection_order = 1, stringsAsFactors = FALSE
  )

  bad_is_qc <- base
  bad_is_qc$is_qc <- "FALSE"
  expect_error(validate_sample_sheet(bad_is_qc), "is_qc.*must be logical")

  bad_injection_order <- base
  bad_injection_order$injection_order <- "1"
  expect_error(validate_sample_sheet(bad_injection_order), "injection_order.*must be numeric")
})

test_that("validate_sample_sheet errors when a referenced file doesn't exist", {
  sheet <- data.frame(
    filepath = "/no/such/file.mzML", column = "RP", polarity = "POS",
    sample_name = "s1", sample_label = "s1", sample_type = "Sample",
    is_qc = FALSE, injection_order = 1, stringsAsFactors = FALSE
  )
  expect_error(validate_sample_sheet(sheet), "don't exist on disk")
})

test_that("validate_sample_sheet passes for a well-formed sheet", {
  f <- tempfile(fileext = ".mzML")
  file.create(f)
  on.exit(unlink(f))

  sheet <- data.frame(
    filepath = f, column = "RP", polarity = "POS", sample_name = "s1",
    sample_label = "s1", sample_type = "Sample", is_qc = FALSE,
    injection_order = 1, stringsAsFactors = FALSE
  )
  expect_true(validate_sample_sheet(sheet))
})
