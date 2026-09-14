test_that("utility functions return expected values", {
  expect_equal(hpsi(c(-2, -1, 0, 2), bend = 1), c(-1, -1, 0, 1))
  expect_equal(keeporder(c("b", "a", "b")), factor(c("b", "a", "b"), levels = c("b", "a")))
  expect_equal(pxly(1:2, 2:3), 0.75)

  paired <- pxlyd(c(1, 3, 2, 5), c(2, 2, 3, 4))
  expect_equal(paired$phat, 0.5)
  expect_length(paired$ci, 2)
  expect_true(all(paired$ci >= 0 & paired$ci <= 1))
})

test_that("quantile summaries respect basic distributional invariants", {
  x <- seq_len(100)

  expect_type(hd(x), "double")
  expect_gt(hd.iqr(x), 0)
  expect_gt(hd.lts(x), 0)
  expect_gt(hd.uts(x), 0)
  expect_equal(hd.crub(x, L = 0, U = 101), (101 - hd(x)) / 101)
  expect_error(hd.crub(x, L = 1), "supplied together")
  expect_error(hd.crub(x, L = 2, U = 1), "L must be")
})

test_that("Winsorized summaries have expected shapes", {
  x <- seq_len(20)
  y <- rev(x)

  expect_type(winvar(x), "double")
  expect_type(trimse(x), "double")
  result <- wincor(x, y)
  expect_named(result, c("n", "cor", "cov", "p.value"))
  expect_equal(result$n, length(x))
})
