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

test_that("qc_tier_available is TRUE when sQC alone has enough", {
  sheet <- make_sheet(2, "sQC")
  expect_true(qc_tier_available(sheet, min_qc = 2))
})

test_that("qc_tier_available is TRUE when only ltQC has enough", {
  sheet <- rbind(
    make_sheet(1, "sQC", filepath = "qc1"),
    make_sheet(2, "ltQC", filepath = paste0("lt", 1:2))
  )
  expect_true(qc_tier_available(sheet, min_qc = 2))
})

test_that("qc_tier_available is FALSE when neither QC tier has enough", {
  sheet <- rbind(
    make_sheet(1, "sQC", filepath = "qc1"),
    make_sheet(1, "ltQC", filepath = "lt1"),
    make_sheet(5, "Sample", filepath = paste0("s", 1:5), injection_order = 1:5)
  )
  expect_false(qc_tier_available(sheet, min_qc = 2))
})

test_that("qc_tier_available excludes qc_flagged rows before counting", {
  sheet <- make_sheet(2, "sQC", qc_flagged = c(FALSE, TRUE))
  expect_false(qc_tier_available(sheet, min_qc = 2))
})

test_that("qc_tier_available treats a missing qc_flagged column as nothing flagged", {
  sheet <- make_sheet(2, "sQC")
  sheet$qc_flagged <- NULL
  expect_true(qc_tier_available(sheet, min_qc = 2))
})

test_that("pick_peaks_cached calls pick_fn and writes a cache on a miss", {
  cache_path <- tempfile()
  on.exit(unlink(cache_path))
  calls <- 0
  stub_pick <- function(filepaths, centwave_param, spectrum_modes) {
    calls <<- calls + 1
    "picked-result"
  }
  result <- pick_peaks_cached(c("a", "b"), list(ppm = 5), c("centroid", "centroid"),
    cache_path = cache_path, pick_fn = stub_pick
  )
  expect_equal(result, "picked-result")
  expect_equal(calls, 1)
  expect_true(file.exists(cache_path))
  expect_false(file.exists(paste0(cache_path, ".tmp")))
})

test_that("pick_peaks_cached reuses the cache when files/params match exactly", {
  cache_path <- tempfile()
  on.exit(unlink(cache_path))
  calls <- 0
  stub_pick <- function(filepaths, centwave_param, spectrum_modes) {
    calls <<- calls + 1
    "picked-result"
  }
  pick_peaks_cached(c("a", "b"), list(ppm = 5), c("centroid", "centroid"),
    cache_path = cache_path, pick_fn = stub_pick
  )
  result <- pick_peaks_cached(c("a", "b"), list(ppm = 5), c("centroid", "centroid"),
    cache_path = cache_path, pick_fn = stub_pick
  )
  expect_equal(result, "picked-result")
  expect_equal(calls, 1) # second call was a cache hit, pick_fn not called again
})

test_that("pick_peaks_cached recomputes when centwave_param differs from the cache", {
  cache_path <- tempfile()
  on.exit(unlink(cache_path))
  calls <- 0
  stub_pick <- function(filepaths, centwave_param, spectrum_modes) {
    calls <<- calls + 1
    sprintf("picked-%d", centwave_param$ppm)
  }
  pick_peaks_cached(c("a", "b"), list(ppm = 5), c("centroid", "centroid"),
    cache_path = cache_path, pick_fn = stub_pick
  )
  result <- pick_peaks_cached(c("a", "b"), list(ppm = 10), c("centroid", "centroid"),
    cache_path = cache_path, pick_fn = stub_pick
  )
  expect_equal(result, "picked-10")
  expect_equal(calls, 2)
})

test_that("pick_peaks_cached recomputes when filepaths differ from the cache", {
  cache_path <- tempfile()
  on.exit(unlink(cache_path))
  calls <- 0
  stub_pick <- function(filepaths, centwave_param, spectrum_modes) {
    calls <<- calls + 1
    "picked-result"
  }
  pick_peaks_cached(c("a", "b"), list(ppm = 5), c("centroid", "centroid"),
    cache_path = cache_path, pick_fn = stub_pick
  )
  pick_peaks_cached(c("a", "c"), list(ppm = 5), c("centroid", "centroid"),
    cache_path = cache_path, pick_fn = stub_pick
  )
  expect_equal(calls, 2)
})

test_that("pick_peaks_cached ignores a corrupt cache file and recomputes", {
  cache_path <- tempfile()
  on.exit(unlink(cache_path))
  writeLines("not a valid rds file", cache_path)
  calls <- 0
  stub_pick <- function(filepaths, centwave_param, spectrum_modes) {
    calls <<- calls + 1
    "picked-result"
  }
  result <- pick_peaks_cached(c("a", "b"), list(ppm = 5), c("centroid", "centroid"),
    cache_path = cache_path, pick_fn = stub_pick
  )
  expect_equal(result, "picked-result")
  expect_equal(calls, 1)
})

test_that("pick_peaks_cached with fresh = TRUE ignores and overwrites an existing cache", {
  cache_path <- tempfile()
  on.exit(unlink(cache_path))
  calls <- 0
  stub_pick <- function(filepaths, centwave_param, spectrum_modes) {
    calls <<- calls + 1
    "picked-result"
  }
  pick_peaks_cached(c("a", "b"), list(ppm = 5), c("centroid", "centroid"),
    cache_path = cache_path, pick_fn = stub_pick
  )
  pick_peaks_cached(c("a", "b"), list(ppm = 5), c("centroid", "centroid"),
    cache_path = cache_path, fresh = TRUE, pick_fn = stub_pick
  )
  expect_equal(calls, 2)
})

test_that("select_center_sample picks the non-flagged sQC closest to the median injection_order", {
  types <- c("sQC", "sQC", "sQC", "ltQC")
  order <- c(1, 5, 10, 3)
  # median(c(1,5,10,3)) == 4 -- order 5 (idx 2) is closer to 4 than 1 or 10.
  result <- select_center_sample(types, order)
  expect_equal(result$idx, 2)
  expect_equal(result$type, "sQC")
})

test_that("select_center_sample falls back to ltQC when no non-flagged sQC is available", {
  types <- c("sQC", "ltQC", "ltQC")
  flagged <- c(TRUE, FALSE, FALSE)
  order <- c(5, 1, 8)
  # median(c(5,1,8)) == 5 -- order 8 (idx 3) is closer to 5 than order 1 is.
  result <- select_center_sample(types, order, qc_flagged = flagged)
  expect_equal(result$idx, 3)
  expect_equal(result$type, "ltQC")
})

test_that("select_center_sample returns NULL when neither tier has any file", {
  types <- c("Sample", "Sample", "Blank")
  order <- c(1, 2, 3)
  result <- select_center_sample(types, order)
  expect_null(result)
})

test_that("select_center_sample excludes qc_flagged rows before picking", {
  types <- c("sQC", "sQC")
  flagged <- c(TRUE, FALSE)
  order <- c(1, 2)
  result <- select_center_sample(types, order, qc_flagged = flagged)
  expect_equal(result$idx, 2)
})

test_that("select_center_sample uses the WHOLE group's median injection_order, not just the candidate pool's", {
  types <- c("sQC", "sQC", "Sample", "Sample", "Sample")
  order <- c(1, 2, 10, 11, 12)
  # Whole-group median(c(1,2,10,11,12)) == 10 -- sQC idx 2 (order 2) is
  # closer to 10 than idx 1 (order 1). A candidate-only median (median of
  # just the two sQC values, 1.5) would wrongly prefer idx 1 instead.
  result <- select_center_sample(types, order)
  expect_equal(result$idx, 2)
})

test_that("select_center_sample treats a missing qc_flagged as nothing flagged", {
  types <- rep("sQC", 2)
  order <- c(1, 2)
  expect_no_error(result <- select_center_sample(types, order))
  expect_true(result$idx %in% c(1, 2))
})

test_that("align_and_correspond errors on an invalid center_sample_mode", {
  skip_if_not_installed("xcms")
  expect_error(
    align_and_correspond(NULL, character(0), center_sample_mode = "bogus"),
    'center_sample_mode must be "qc" or "middle"'
  )
})
