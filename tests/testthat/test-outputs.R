test_that("parametric inference returns stable result structures", {
  x <- seq_len(30)
  y <- x + rep(c(-1, 1), length.out = length(x))

  trim_result <- trimci(x)
  yuen_result <- yuen(x, y)
  yuend_result <- yuend(x, y)

  expect_named(trim_result, c("ci", "estimate", "test.stat", "se", "df", "p.value", "alternative", "n"))
  expect_named(yuen_result, c("n1", "n2", "est.1", "est.2", "ci", "p.value", "dif", "se", "teststat", "crit", "df", "alternative"))
  expect_named(yuend_result, c("ci", "p.value", "est1", "est2", "dif", "se", "teststat", "n", "df"))

  expect_length(trim_result$ci, 2)
  expect_length(yuen_result$ci, 2)
  expect_length(yuend_result$ci, 2)
  expect_true(all(trim_result$ci == sort(trim_result$ci)))
  expect_true(all(yuen_result$ci == sort(yuen_result$ci)))
  expect_true(all(yuend_result$ci == sort(yuend_result$ci)))
  expect_true(all(c(trim_result$p.value, yuen_result$p.value, yuend_result$p.value) >= 0))
  expect_true(all(c(trim_result$p.value, yuen_result$p.value, yuend_result$p.value) <= 1))
})

test_that("percentile bootstrap output retains inference metadata", {
  set.seed(42)
  result <- twosampb(rnorm(20), rnorm(25), nboot = 40)

  expect_named(result, c("est.1", "est.2", "estimate", "ci", "p.value", "hyp", "alternative", "small.n", "sq.se", "boot.estimates", "boot.estimates.x", "boot.estimates.y", "ind", "n1", "n2"))
  expect_length(result$boot.estimates, 40)
  expect_equal(result$n1, 20)
  expect_equal(result$n2, 25)
  expect_true(isTRUE(result$ind))
})
