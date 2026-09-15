test_that("plot methods return ggplot objects for supported inputs", {
  set.seed(5)

  expect_s3_class(plot_boot(rnorm(100)), "ggplot")
  expect_s3_class(
    plot_boot(list(boot.estimates = rnorm(100), estimate = 0)),
    "ggplot"
  )
  expect_s3_class(plot_boot2(rnorm(100), rnorm(100)), "ggplot")
  expect_s3_class(
    plot_boot2(list(boot.estimates.x = rnorm(100), boot.estimates.y = rnorm(100))),
    "ggplot"
  )
  expect_s3_class(plot_ecdf(rnorm(20), rnorm(20)), "ggplot")
})

test_that("plot methods validate required numeric inputs", {
  expect_error(plot_boot(1), "at least two")
  expect_error(plot_boot2(1:3, 1:2), "equal length")
})

test_that("plot layers accept additional ggplot2 arguments", {
  expect_s3_class(plot_boot2(rnorm(100), rnorm(100), shape = 21, size = 2), "ggplot")
  expect_s3_class(plot_ecdf(rnorm(20), rnorm(20), na.rm = TRUE), "ggplot")
})
