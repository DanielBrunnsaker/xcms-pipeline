test_that("report_qc_batch_coverage is silent when every batch has both QC types", {
  sheet <- data.frame(
    batch = c("B1", "B1", "B2", "B2"),
    sample_type = c("sQC", "ltQC", "sQC", "ltQC"),
    qc_flagged = FALSE,
    stringsAsFactors = FALSE
  )
  expect_no_warning(report_qc_batch_coverage(sheet))
})

test_that("report_qc_batch_coverage warns when a batch has no sQC at all", {
  sheet <- data.frame(
    batch = c("B1", "B1", "B2", "B2"),
    sample_type = c("sQC", "ltQC", "ltQC", "ltQC"),
    qc_flagged = FALSE,
    stringsAsFactors = FALSE
  )
  expect_warning(report_qc_batch_coverage(sheet), "sQC has no non-flagged files.*B2")
})

test_that("report_qc_batch_coverage treats a fully-flagged batch as missing coverage", {
  sheet <- data.frame(
    batch = c("B1", "B1", "B2"),
    sample_type = c("sQC", "ltQC", "sQC"),
    qc_flagged = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  expect_warning(report_qc_batch_coverage(sheet), "sQC has no non-flagged files.*B2")
})

test_that("report_qc_batch_coverage treats a missing qc_flagged column as nothing flagged", {
  sheet <- data.frame(
    batch = c("B1", "B1", "B2", "B2"),
    sample_type = c("sQC", "ltQC", "sQC", "ltQC"),
    stringsAsFactors = FALSE
  )
  expect_no_warning(report_qc_batch_coverage(sheet))
})

test_that("report_qc_batch_coverage warns separately for both QC types when both are incomplete", {
  sheet <- data.frame(
    batch = c("B1", "B2"),
    sample_type = c("sQC", "ltQC"),
    qc_flagged = FALSE,
    stringsAsFactors = FALSE
  )
  warnings <- character(0)
  withCallingHandlers(
    report_qc_batch_coverage(sheet),
    warning = function(w) {
      warnings[[length(warnings) + 1]] <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("sQC has no non-flagged files.*B2", warnings)))
  expect_true(any(grepl("ltQC has no non-flagged files.*B1", warnings)))
})
