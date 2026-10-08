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

test_that("joint two-sample bootstrap inference returns expected outputs", {
  set.seed(5)
  x <- rnorm(20, mean = 1)
  y <- rnorm(25)

  result <- twosampb_joint(x, y, nboot = 40)

  expect_equal(result$estimate, mxmy(x, y))
  expect_length(result$boot.estimates, 40)
  expect_equal(result$n1, length(x))
  expect_equal(result$n2, length(y))
  expect_true(all(is.finite(result$ci)))
  expect_true(result$p.value >= 0 && result$p.value <= 1)

  result_pxly <- twosampb_joint(x, y, est = pxly, hyp = 0.5, nboot = 40)
  expect_equal(result_pxly$estimate, pxly(x, y))
})

test_that("joint two-sample bootstrap supports paired data and one-sided intervals", {
  set.seed(6)
  x <- rnorm(20)
  y <- x + rnorm(20)

  paired <- twosampb_joint(cbind(x, y), ind = FALSE, nboot = 40)
  expect_false(paired$ind)
  expect_equal(paired$n1, paired$n2)

  greater <- twosampb_joint(x, y, alternative = "greater", nboot = 40)
  expect_identical(greater$ci[2], Inf)
  less <- twosampb_joint(x, y, alternative = "less", nboot = 40)
  expect_identical(less$ci[1], -Inf)
})

test_that("joint two-sample bootstrap validates inputs", {
  x <- rnorm(10)
  expect_error(twosampb_joint(x, nboot = 10), "y must be supplied")
  expect_error(twosampb_joint(x, x[-1], ind = FALSE, nboot = 10), "same length")
  expect_error(twosampb_joint(x, x, small.n = NA, nboot = 10), "small.n")
})

test_that("Yuenbt bootstrap inference returns expected outputs for one-sided confidence intervals", {
  set.seed(4)
  result <- yuenbt(rnorm(25, mean = 1), rnorm(25),
                   nboot = 40, alternative = "greater")

  expect_identical(result$alternative, "greater")
  expect_true(is.finite(result$ci[1]))
  expect_identical(result$ci[2], Inf)
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})
