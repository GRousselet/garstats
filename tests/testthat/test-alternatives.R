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

expect_wider_ci <- function(narrower, wider) {
  expect_lte(wider[1], narrower[1])
  expect_gte(wider[2], narrower[2])
}

expect_one_sided_ci <- function(result, alternative) {
  if (alternative == "greater") {
    expect_true(is.finite(result$ci[1]))
    expect_identical(result$ci[2], Inf)
  } else {
    expect_identical(result$ci[1], -Inf)
    expect_true(is.finite(result$ci[2]))
  }
  expect_true(is.finite(result$p.value))
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
}

test_that("confidence intervals widen as alpha decreases", {
  x <- seq(-5, 24)
  y <- seq(-7, 22)

  expect_wider_ci(
    trimci(x, alpha = 0.05)$ci,
    trimci(x, alpha = 0.01)$ci
  )
  expect_wider_ci(
    yuen(x, y, alpha = 0.05)$ci,
    yuen(x, y, alpha = 0.01)$ci
  )
  expect_wider_ci(
    yuend(x, y, alpha = 0.05)$ci,
    yuend(x, y, alpha = 0.01)$ci
  )
  expect_wider_ci(
    sint(x, alpha = 0.05),
    sint(x, alpha = 0.01)
  )
  expect_wider_ci(
    pxlyd(x, y, alpha = 0.05)$ci,
    pxlyd(x, y, alpha = 0.01)$ci
  )
})

test_that("bootstrap confidence intervals widen as alpha decreases", {
  x <- seq(-5, 24)
  y <- seq(-7, 22)
  set.seed(10)
  bootstrap_x <- rnorm(30)
  bootstrap_y <- bootstrap_x + rnorm(30, sd = 0.5)

  set.seed(11)
  trimcibt_05 <- trimcibt(bootstrap_x, alpha = 0.05, nboot = 80)$ci
  set.seed(11)
  trimcibt_01 <- trimcibt(bootstrap_x, alpha = 0.01, nboot = 80)$ci
  expect_wider_ci(trimcibt_05, trimcibt_01)

  set.seed(12)
  yuenbt_05 <- yuenbt(bootstrap_x, bootstrap_y, alpha = 0.05, nboot = 80)$ci
  set.seed(12)
  yuenbt_01 <- yuenbt(bootstrap_x, bootstrap_y, alpha = 0.01, nboot = 80)$ci
  expect_wider_ci(yuenbt_05, yuenbt_01)

  set.seed(13)
  yuendbt_05 <- yuendbt(bootstrap_x, bootstrap_y, alpha = 0.05, nboot = 80)$ci
  set.seed(13)
  yuendbt_01 <- yuendbt(bootstrap_x, bootstrap_y, alpha = 0.01, nboot = 80)$ci
  expect_wider_ci(yuendbt_05, yuendbt_01)

  set.seed(14)
  onesampb_05 <- onesampb(bootstrap_x, alpha = 0.05, nboot = 80)$ci
  set.seed(14)
  onesampb_01 <- onesampb(bootstrap_x, alpha = 0.01, nboot = 80)$ci
  expect_wider_ci(onesampb_05, onesampb_01)

  set.seed(15)
  twosampb_05 <- twosampb(bootstrap_x, bootstrap_y, alpha = 0.05, nboot = 80)$ci
  set.seed(15)
  twosampb_01 <- twosampb(bootstrap_x, bootstrap_y, alpha = 0.01, nboot = 80)$ci
  expect_wider_ci(twosampb_05, twosampb_01)

  set.seed(16)
  bootstrap_estimates <- rnorm(80)
  expect_wider_ci(
    pbci(bootstrap_estimates, alpha = 0.05)$ci,
    pbci(bootstrap_estimates, alpha = 0.01)$ci
  )
})

test_that("one-sided confidence intervals use the requested alternative", {
  x <- seq(-5, 24)
  y <- seq(-7, 22)
  set.seed(17)
  bootstrap_x <- rnorm(30)
  bootstrap_y <- bootstrap_x + rnorm(30, sd = 0.5)

  for (alternative in c("greater", "less")) {
    expect_one_sided_ci(trimci(x, alternative = alternative), alternative)
    expect_one_sided_ci(
      trimcibt(bootstrap_x, alternative = alternative, nboot = 80),
      alternative
    )
    expect_one_sided_ci(yuen(x, y, alternative = alternative), alternative)
    expect_one_sided_ci(
      yuenbt(bootstrap_x, bootstrap_y, alternative = alternative, nboot = 80),
      alternative
    )
    expect_one_sided_ci(
      yuendbt(bootstrap_x, bootstrap_y, alternative = alternative, nboot = 80),
      alternative
    )
    expect_one_sided_ci(
      onesampb(bootstrap_x, alternative = alternative, nboot = 80),
      alternative
    )
    expect_one_sided_ci(
      twosampb(bootstrap_x, bootstrap_y, alternative = alternative, nboot = 80),
      alternative
    )
    expect_one_sided_ci(
      pbci(rnorm(80), alternative = alternative),
      alternative
    )
  }
})
