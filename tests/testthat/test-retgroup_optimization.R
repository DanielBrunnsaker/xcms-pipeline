test_that("report_qc_batch_coverage prints no WARNING when every batch has both QC types", {
  sheet <- data.frame(
    batch = c("B1", "B1", "B2", "B2"),
    sample_type = c("sQC", "ltQC", "sQC", "ltQC"),
    qc_flagged = FALSE,
    stringsAsFactors = FALSE
  )
  msgs <- testthat::capture_messages(report_qc_batch_coverage(sheet))
  expect_false(any(grepl("WARNING:", msgs, fixed = TRUE)))
})

test_that("report_qc_batch_coverage flags a batch with no sQC at all", {
  sheet <- data.frame(
    batch = c("B1", "B1", "B2", "B2"),
    sample_type = c("sQC", "ltQC", "ltQC", "ltQC"),
    qc_flagged = FALSE,
    stringsAsFactors = FALSE
  )
  msgs <- testthat::capture_messages(report_qc_batch_coverage(sheet))
  expect_true(any(grepl("WARNING: sQC has no non-flagged files.*B2", msgs)))
})

test_that("report_qc_batch_coverage treats a fully-flagged batch as missing coverage", {
  sheet <- data.frame(
    batch = c("B1", "B1", "B2"),
    sample_type = c("sQC", "ltQC", "sQC"),
    qc_flagged = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  msgs <- testthat::capture_messages(report_qc_batch_coverage(sheet))
  expect_true(any(grepl("WARNING: sQC has no non-flagged files.*B2", msgs)))
})

test_that("report_qc_batch_coverage treats a missing qc_flagged column as nothing flagged", {
  sheet <- data.frame(
    batch = c("B1", "B1", "B2", "B2"),
    sample_type = c("sQC", "ltQC", "sQC", "ltQC"),
    stringsAsFactors = FALSE
  )
  msgs <- testthat::capture_messages(report_qc_batch_coverage(sheet))
  expect_false(any(grepl("WARNING:", msgs, fixed = TRUE)))
})

test_that("report_qc_batch_coverage flags both QC types separately when both are incomplete", {
  sheet <- data.frame(
    batch = c("B1", "B2"),
    sample_type = c("sQC", "ltQC"),
    qc_flagged = FALSE,
    stringsAsFactors = FALSE
  )
  msgs <- testthat::capture_messages(report_qc_batch_coverage(sheet))
  expect_true(any(grepl("WARNING: sQC has no non-flagged files.*B2", msgs)))
  expect_true(any(grepl("WARNING: ltQC has no non-flagged files.*B1", msgs)))
})
