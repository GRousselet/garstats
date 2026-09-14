test_that("one-sample bootstrap inference returns expected outputs", {
  set.seed(2)
  x <- rnorm(25)

  result <- trimcibt(x, nboot = 40)

  expect_length(result$boot.estimates, 40)
  expect_length(result$boot.t, 40)
  expect_equal(result$n, length(x))
  expect_true(all(is.finite(result$ci)))
})

test_that("two-sample bootstrap inference returns expected outputs", {
  set.seed(3)
  x <- rnorm(25)
  y <- x + rnorm(25, sd = 0.5)

  result <- twosampb(x, y, ind = FALSE, nboot = 40)

  expect_false(result$ind)
  expect_equal(result$n1, length(x))
  expect_equal(result$n2, length(y))
  expect_length(result$boot.estimates, 40)
  expect_length(result$boot.estimates.x, 40)
  expect_length(result$boot.estimates.y, 40)
})

test_that("bootstrap results honor one-sided confidence intervals", {
  set.seed(4)
  result <- yuenbt(rnorm(25, mean = 1), rnorm(25),
                   nboot = 40, alternative = "greater")

  expect_identical(result$alternative, "greater")
  expect_true(is.finite(result$ci[1]))
  expect_identical(result$ci[2], Inf)
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})
