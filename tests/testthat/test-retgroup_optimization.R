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

test_that("select_retgroup_qc_idx picks sQC automatically when there are enough", {
  types <- c("sQC", "sQC", "ltQC")
  result <- select_retgroup_qc_idx(types, min_qc = 2)
  expect_equal(result$type, "sQC")
  expect_false(result$forced)
  expect_equal(result$idx, c(1, 2))
})

test_that("select_retgroup_qc_idx falls back to ltQC automatically when sQC is insufficient", {
  types <- c("sQC", "ltQC", "ltQC")
  result <- select_retgroup_qc_idx(types, min_qc = 2)
  expect_equal(result$type, "ltQC")
  expect_false(result$forced)
  expect_equal(result$idx, c(2, 3))
})

test_that("select_retgroup_qc_idx errors when neither type has enough", {
  types <- c("sQC", "ltQC", "Sample")
  expect_error(select_retgroup_qc_idx(types, min_qc = 2), "Not enough non-flagged QC files")
})

test_that("select_retgroup_qc_idx excludes qc_flagged rows before counting", {
  types <- c("sQC", "sQC", "ltQC", "ltQC")
  flagged <- c(TRUE, FALSE, FALSE, FALSE)
  result <- select_retgroup_qc_idx(types, qc_flagged = flagged, min_qc = 2)
  # Only 1 non-flagged sQC left (below min_qc) -> falls back to ltQC.
  expect_equal(result$type, "ltQC")
  expect_equal(result$idx, c(3, 4))
})

test_that("select_retgroup_qc_idx honors a forced qc_type even when the automatic pick would differ", {
  types <- c("sQC", "sQC", "sQC", "ltQC", "ltQC")
  result <- select_retgroup_qc_idx(types, qc_type = "ltQC", min_qc = 2)
  expect_equal(result$type, "ltQC")
  expect_true(result$forced)
  expect_equal(result$idx, c(4, 5))
})

test_that("select_retgroup_qc_idx forced qc_type is case-insensitive", {
  types <- c("sQC", "sQC", "ltQC", "ltQC")
  result <- select_retgroup_qc_idx(types, qc_type = "SQC", min_qc = 2)
  expect_equal(result$type, "sQC")
  expect_true(result$forced)
})

test_that("select_retgroup_qc_idx errors if the forced qc_type doesn't have enough files", {
  types <- c("sQC", "sQC", "sQC", "ltQC")
  expect_error(
    select_retgroup_qc_idx(types, qc_type = "ltQC", min_qc = 2),
    "qc_type forced to ltQC, but only 1"
  )
})

test_that("select_retgroup_qc_idx errors on an unrecognized forced qc_type", {
  types <- c("sQC", "sQC")
  expect_error(select_retgroup_qc_idx(types, qc_type = "blank", min_qc = 2), 'must be "sQC" or "ltQC"')
})
