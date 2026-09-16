test_that("trimmed-mean functions validate trimming proportions", {
  x <- seq_len(20)

  expect_error(trimci(x, tr = -0.1), "tr must be")
  expect_error(trimse(x, tr = 0.5), "tr must be")
  expect_error(yuen(x, x + 1, tr = 0.5), "tr must be")
  expect_error(yuenbt(x, x + 1, tr = Inf, nboot = 10), "tr must be")
})

test_that("paired functions require equal sample sizes", {
  x <- seq_len(20)

  expect_error(sint(x, x[-1]), "equal sample sizes")
  expect_error(yuend(x, x[-1]), "equal sample sizes")
  expect_error(yuendbt(x, x[-1], nboot = 10), "equal sample sizes")
  expect_error(
    twosampb(x, x[-1], ind = FALSE, nboot = 10),
    "same length"
  )
})

test_that("bootstrap interfaces reject invalid options", {
  x <- seq_len(20)

  expect_error(sint(x, alpha = 0), "alpha must be")
  expect_error(sint(x, alpha = 1), "alpha must be")
  expect_error(trimcibt(x, small.n = NA, nboot = 10), "small.n")
  expect_error(pbci("not numeric"), "numeric vector")
  expect_error(
    twosampb(x, ind = FALSE, nboot = 10),
    "matrix"
  )
})
