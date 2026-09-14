test_that("Yuen alternatives have matching p-values and confidence intervals", {
  set.seed(1)
  x <- rnorm(30, mean = 1)
  y <- rnorm(35)

  two_sided <- yuen(x, y, alternative = "two.sided")
  greater <- yuen(x, y, alternative = "greater")
  less <- yuen(x, y, alternative = "less")

  expect_true(all(is.finite(two_sided$ci)))
  expect_true(is.finite(greater$ci[1]))
  expect_identical(greater$ci[2], Inf)
  expect_identical(less$ci[1], -Inf)
  expect_true(is.finite(less$ci[2]))
  expect_equal(
    two_sided$p.value,
    min(1, 2 * min(greater$p.value, less$p.value)),
    tolerance = 1e-12
  )
  expect_lt(greater$p.value, less$p.value)
})

test_that("smaller alpha yields more conservative Yuen confidence bounds", {
  x <- seq(-5, 24)
  y <- seq(-7, 22)

  two_95 <- yuen(x, y, alpha = 0.05)$ci
  two_99 <- yuen(x, y, alpha = 0.01)$ci
  greater_95 <- yuen(x, y, alpha = 0.05, alternative = "greater")$ci
  greater_99 <- yuen(x, y, alpha = 0.01, alternative = "greater")$ci
  less_95 <- yuen(x, y, alpha = 0.05, alternative = "less")$ci
  less_99 <- yuen(x, y, alpha = 0.01, alternative = "less")$ci

  expect_lte(two_99[1], two_95[1])
  expect_gte(two_99[2], two_95[2])
  expect_lte(greater_99[1], greater_95[1])
  expect_gte(less_99[2], less_95[2])
})

test_that("trimci alternatives follow the same confidence-bound convention", {
  x <- seq(-10, 19)

  expect_identical(trimci(x, alternative = "greater")$ci[2], Inf)
  expect_identical(trimci(x, alternative = "less")$ci[1], -Inf)
})
