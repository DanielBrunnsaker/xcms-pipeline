test_that("low_outlier_threshold is a fraction of the median", {
  x <- c(10, 20, 30, 40, 50)
  expect_equal(low_outlier_threshold(x, 0.5), 15)
  expect_equal(low_outlier_threshold(x, 1), 30)
})

test_that("to_rgba formats a named color with the given alpha", {
  expect_equal(to_rgba("black", 1), "rgba(0,0,0,1.00)")
  expect_equal(to_rgba("white", 0.3), "rgba(255,255,255,0.30)")
})

test_that("compute_type_thresholds computes each type's threshold from only its own values", {
  # sQC much lower than ltQC -- e.g. ltQC more concentrated, as observed in
  # practice. Each type's threshold should reflect only its own median, not
  # a pooled one.
  x <- c(100, 120, 80, 500, 520, 480)
  types <- c("sQC", "sQC", "sQC", "ltQC", "ltQC", "ltQC")
  result <- compute_type_thresholds(x, types, usable_idx = 1:6, min_fraction_of_median = 0.5)
  expect_equal(result[["sQC"]], 50)   # 0.5 * median(100,120,80) = 0.5*100
  expect_equal(result[["ltQC"]], 250) # 0.5 * median(500,520,480) = 0.5*500
})

test_that("compute_type_thresholds restricts to usable_idx", {
  x <- c(100, 5, 120, 80)
  types <- c("sQC", "sQC", "sQC", "sQC")
  # Row 2 (value 5) excluded from usable_idx -- shouldn't drag the median down.
  result <- compute_type_thresholds(x, types, usable_idx = c(1, 3, 4), min_fraction_of_median = 0.5)
  expect_equal(result[["sQC"]], 50) # 0.5 * median(100, 120, 80) = 0.5*100
})

test_that("compute_type_thresholds gives NA for a type with zero usable rows", {
  x <- c(100, 120, 500)
  types <- c("sQC", "sQC", "ltQC")
  # ltQC's only row (index 3) excluded from usable_idx.
  result <- compute_type_thresholds(x, types, usable_idx = c(1, 2), min_fraction_of_median = 0.5)
  expect_equal(result[["sQC"]], 55) # 0.5 * median(100, 120) = 0.5*110
  expect_true(is.na(result[["ltQC"]]))
})

test_that("compute_type_thresholds handles a single sample_type (e.g. sQC-only checks)", {
  x <- c(10, 20, 30, 40, 50)
  types <- rep("sQC", 5)
  result <- compute_type_thresholds(x, types, usable_idx = 1:5, min_fraction_of_median = 0.5)
  expect_equal(unname(result), 15)
  expect_equal(names(result), "sQC")
})

test_that("a lower min_fraction_of_median gives a more lenient (lower) threshold", {
  # Simulates check_qc_quality()'s much-more-lenient
  # aligned_min_fraction_of_median vs. the stricter raw-peak-count one.
  x <- c(80, 90, 100, 110, 120)
  types <- rep("sQC", 5)
  strict <- compute_type_thresholds(x, types, usable_idx = 1:5, min_fraction_of_median = 0.5)
  lenient <- compute_type_thresholds(x, types, usable_idx = 1:5, min_fraction_of_median = 0.2)
  expect_lt(lenient[["sQC"]], strict[["sQC"]])
  # A file at 30% of the median falls below the strict (50%) threshold but
  # not the lenient (20%) one -- the exact scenario a technically-fine-but-
  # poorly-aligned QC should trigger only the secondary (feature-count)
  # check, not both.
  value <- 0.3 * median(x)
  expect_true(value < strict[["sQC"]])
  expect_false(value < lenient[["sQC"]])
})
