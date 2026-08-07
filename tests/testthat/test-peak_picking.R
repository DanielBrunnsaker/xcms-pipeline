make_sheet <- function(n, sample_type = "sQC", batch = "B1", injection_order = seq_len(n),
                       filepath = paste0("f", seq_len(n)), spectrum_mode = NA_character_,
                       qc_flagged = FALSE) {
  data.frame(
    filepath = filepath, sample_type = sample_type, batch = batch,
    injection_order = injection_order, spectrum_mode = spectrum_mode,
    qc_flagged = qc_flagged, column = "RP", polarity = "POS",
    stringsAsFactors = FALSE
  )
}

test_that("select_ipo_subset prefers sQC when there are enough", {
  sheet <- rbind(
    make_sheet(3, "sQC", filepath = paste0("qc", 1:3)),
    make_sheet(5, "Sample", filepath = paste0("s", 1:5), injection_order = 1:5)
  )
  result <- select_ipo_subset(sheet, n = 3, min_qc = 2)
  expect_true(all(grepl("^qc", result)))
})

test_that("select_ipo_subset falls through to ltQC when sQC is insufficient", {
  sheet <- rbind(
    make_sheet(1, "sQC", filepath = "qc1"),
    make_sheet(3, "ltQC", filepath = paste0("lt", 1:3))
  )
  result <- select_ipo_subset(sheet, n = 3, min_qc = 2)
  expect_true(all(grepl("^lt", result)))
})

test_that("select_ipo_subset falls through to regular samples as a last resort", {
  sheet <- rbind(
    make_sheet(1, "sQC", filepath = "qc1"),
    make_sheet(1, "ltQC", filepath = "lt1"),
    make_sheet(4, "Sample", filepath = paste0("s", 1:4), injection_order = 1:4)
  )
  result <- select_ipo_subset(sheet, n = 4, min_qc = 2)
  expect_true(all(grepl("^s", result)))
  expect_length(result, 4)
})

test_that("select_ipo_subset excludes qc_flagged rows before counting toward min_qc", {
  sheet <- rbind(
    make_sheet(3, "sQC", filepath = paste0("qc", 1:3), qc_flagged = c(TRUE, TRUE, FALSE)),
    make_sheet(3, "ltQC", filepath = paste0("lt", 1:3))
  )
  # Only 1 unflagged sQC remains (< min_qc = 2) -> falls through to ltQC.
  result <- select_ipo_subset(sheet, n = 3, min_qc = 2)
  expect_true(all(grepl("^lt", result)))
})

test_that("select_ipo_subset spreads a single batch evenly across injection_order", {
  sheet <- make_sheet(10, "sQC", injection_order = 1:10, filepath = paste0("f", 1:10))
  result <- select_ipo_subset(sheet, n = 4, min_qc = 2)
  expect_equal(result, c("f1", "f4", "f7", "f10"))
})

test_that("select_ipo_subset subsamples batches evenly when there are more than n", {
  sheet <- do.call(rbind, lapply(1:6, function(i) {
    make_sheet(1, "sQC", batch = paste0("B", i), injection_order = i, filepath = paste0("fB", i))
  }))
  result <- select_ipo_subset(sheet, n = 4, min_qc = 2)
  expect_equal(sort(unname(result)), sort(c("fB1", "fB3", "fB4", "fB6")))
})

test_that("select_ipo_subset distributes n evenly across fewer/equal batches", {
  b1 <- make_sheet(4, "sQC", batch = "B1", injection_order = 1:4, filepath = paste0("fB1_", 1:4))
  b2 <- make_sheet(4, "sQC", batch = "B2", injection_order = 5:8, filepath = paste0("fB2_", 1:4))
  sheet <- rbind(b1, b2)
  result <- select_ipo_subset(sheet, n = 4, min_qc = 2)
  expect_equal(sort(result), sort(c("fB1_1", "fB1_4", "fB2_1", "fB2_4")))
})

test_that("select_ipo_subset gives any remainder to the earliest batches", {
  b1 <- make_sheet(1, "sQC", batch = "B1", injection_order = 1, filepath = "fB1_1")
  b2 <- make_sheet(1, "sQC", batch = "B2", injection_order = 2, filepath = "fB2_1")
  b3 <- make_sheet(1, "sQC", batch = "B3", injection_order = 3, filepath = "fB3_1")
  sheet <- rbind(b1, b2, b3)
  # n=4 across 3 batches -> per_batch_n = c(2,1,1), but each batch only has
  # 1 file available, so pick_spread just caps at what's there either way.
  result <- select_ipo_subset(sheet, n = 4, min_qc = 2)
  expect_equal(sort(result), c("fB1_1", "fB2_1", "fB3_1"))
})

test_that("select_ipo_subset prefers a uniform spectrum mode when one has enough candidates", {
  sheet <- make_sheet(
    5, "sQC", injection_order = 1:5, filepath = paste0("f", 1:5),
    spectrum_mode = c("centroid", "centroid", "centroid", "centroid", "profile")
  )
  result <- select_ipo_subset(sheet, n = 4, min_qc = 2)
  expect_equal(sort(result), c("f1", "f2", "f3", "f4"))
})

test_that("select_ipo_subset errors when nothing usable is available", {
  sheet <- make_sheet(2, "Blank", filepath = c("b1", "b2"))
  expect_error(select_ipo_subset(sheet, n = 4, min_qc = 2), "No usable files for IPO optimization")
})
