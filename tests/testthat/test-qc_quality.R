test_that("low_outlier_threshold is a fraction of the median", {
  x <- c(10, 20, 30, 40, 50)
  expect_equal(low_outlier_threshold(x, 0.5), 15)
  expect_equal(low_outlier_threshold(x, 1), 30)
})

test_that("flag_low_outliers only flags values below the threshold", {
  x <- c(5, 10, 20, 30, 40)
  # median = 20, threshold at 0.5 = 10 -> only 5 is strictly below it.
  expect_equal(flag_low_outliers(x, 0.5), c(TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("flag_low_outliers never flags high values", {
  # median = 50, threshold = 25 -- neither the typical values nor the high
  # outlier fall below it.
  x <- c(50, 50, 50, 1000)
  expect_equal(flag_low_outliers(x, 0.5), c(FALSE, FALSE, FALSE, FALSE))
})

test_that("to_rgba formats a named color with the given alpha", {
  expect_equal(to_rgba("black", 1), "rgba(0,0,0,1.00)")
  expect_equal(to_rgba("white", 0.3), "rgba(255,255,255,0.30)")
})
