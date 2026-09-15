# ===================================================
# The functions in this file were modified from the
# original code from Rand R. Wilcox:
# http://dornsife.usc.edu/labs/rwilcox/software/
# https://osf.io/xhe8u/
# The main reference for the code is:
# Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
# Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
# ===================================================

# MAJOR CHANGES
# - remove calls to elimna(); it is better practice to remove missing values before calling a function
# x <- na.omit(x) # vector
# x <- as.vector(na.omit(x)) # to get a plain vector
# m <- m[complete.cases(m), , drop = FALSE] # matrix or dataframe
# x <- lapply(x, na.omit) # list of vectors
# - change default bootstrap samples to 2000
# - give user the option to return two-sided or one-sided p-values and confidence intervals
# - remove option to call set.seed() in the functions; it is better practice to set the seed outside of the functions
# - use quantile() to compute quantiles -- default type=6
# - export bootstrap distributions
# - boostrap p-values are by default corrected for small sample sizes

#' @importFrom stats cor mad median optim pbeta pbinom pt qbinom qnorm qt quantile var
NULL

#' Harrell-Davis quantile estimator
#'
#' @param x Numeric data vector.
#' @param q Desired quantile, between 0 and 1.
#' @return A numeric Harrell-Davis quantile estimate.
#' @export
#'
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
#'
#' Harrell, F. E., & Davis, C. E. (1982). A new distribution-free quantile estimator.
#' Biometrika, 69(3), 635–640. https://doi.org/10.1093/biomet/69.3.635
hd <- function(x, q = 0.5) {
  n <- length(x)
  m1 <- (n + 1) * q
  m2 <- (n + 1) * (1 - q)
  vec <- seq(along = x)
  w <- pbeta(vec / n, m1, m2) - pbeta((vec - 1) / n, m1, m2) # W sub i values
  y <- sort(x)
  hd <- sum(w * y)
  hd
}

#' Parametric inference for a one-sample trimmed mean
#'
#' @param x Numeric data vector. Missing values must be removed first.
#' @param tr Trimming proportion in [0, 0.5), with 0.2 default.
#' @param alpha Significance level, with 0.03 default (97% confidence).
#' @param hyp (Null) hypothesis value, with 0 default.
#' @param alternative Alternative hypothesis: "two.sided", "greater", or
#'   "less", matching the alternative argument in t.test().
#'
#' @return A list containing the confidence interval, estimate, test statistic,
#'   standard error, degrees of freedom, and p-value.
#'
#' @export
#'
#' @details
#' For "greater" the confidence interval is `[lower bound, Inf]`
#' For "less" the confidence interval is `[-Inf, upper bound]`.
#'
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
#'
#' @examples
#' n <- 20
#' set.seed(666)
#' # No effect: population trimmed mean = 0
#' x <- rnorm(20, mean = 0)
#' trimci(x, alternative = "two.sided")
#' trimci(x, alternative = "greater")
#' trimci(x, alternative = "less")
#' # Negative effect: population trimmed mean = -1
#' x <- rnorm(20, mean = -1)
#' trimci(x, alternative = "less")
#' # Positive effect: population trimmed mean = 1
#' x <- rnorm(20, mean = 1)
#' trimci(x, alternative = "greater")
#' # Inferences on means and 20% trimmed means return very different results:
#' x <- rlnorm(50)
#' trimci(x, tr = 0, alternative = "less", hyp = 1.5)
#' trimci(x, tr = 0.2, alternative = "less", hyp = 1.5)
trimci <- function(x, tr = .2, alpha = .03, hyp = 0,
                   alternative = c("two.sided", "greater", "less")) {
  if (!is.numeric(tr) || length(tr) != 1L || !is.finite(tr) || tr < 0 || tr >= 0.5) {
    stop("tr must be a single finite number in [0, 0.5).")
  }
  alternative <- match.arg(alternative)
  n <- length(x)
  estimate <- mean(x, tr)
  se <- sqrt(winvar(x, tr)) /
    ((1 - 2 * tr) * sqrt(n))
  df <- n - 2 * floor(tr * n) - 1
  test.stat <- (estimate - hyp) / se
  if (alternative == "two.sided") { # two-sided confidence interval and p-value
    critical <- qt(1 - alpha / 2, df)
    ci <- estimate + c(-1, 1) * critical * se
    p.value <- 2 * pt(-abs(test.stat), df)
  } else { # one-sided confidence interval and p-value
    critical <- qt(1 - alpha, df)
    if (alternative == "greater") {
      ci <- c(estimate - critical * se, Inf)
      p.value <- pt(test.stat, df, lower.tail = FALSE)
    } else if (alternative == "less") {
      ci <- c(-Inf, estimate + critical * se)
      p.value <- pt(test.stat, df, lower.tail = TRUE)
    }
  }
  list(
    ci = ci,
    estimate = estimate,
    test.stat = test.stat,
    se = se,
    df = df,
    p.value = p.value,
    alternative = alternative,
    n = n
  )
}

#' One-sample bootstrap-t inference for a trimmed mean
#'
#' Computes a bootstrap percentile-t confidence interval for a trimmed mean.
#' Missing values are not removed internally. Remove them before calling this
#' function.
#'
#' @param x Numeric data vector.
#' @param tr Trimming proportion in [0, 0.5).
#' @param alpha Significance level.
#' @param nboot Number of bootstrap samples.
#' @param ci.type Type of confidence interval for two-sided inference: "symmetric" or
#'   "asymmetric". Ignored for one-sided inference, which always uses signed
#'   bootstrap statistics.
#' @param q.type Type of quantile estimator to compute bootstrap interval. Default to 6.
#' @param hyp Hypothesis value for the trimmed mean. Default to the null.
#' @param alternative Alternative hypothesis: "two.sided", "greater", or
#'   "less".
#' @param small.n Logical; apply the small-sample correction
#'   `(tail.count + 1) / (nboot + 1)` to the bootstrap p-value.
#'
#' @return A list containing the trimmed-mean estimate, confidence interval,
#'   test statistic, p-value, sample size, bootstrap estimates, and bootstrap
#'   t-statistics.
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
#'
#' @examples
#' set.seed(666)
#' x <- rnorm(20, mean = 1)
#' trimcibt(x, ci.type = "symmetric")
#' trimcibt(x, ci.type = "asymmetric")
#' trimcibt(x, ci.type = "asymmetric", alternative = "greater")
#' trimcibt(x, ci.type = "asymmetric", alternative = "less")
#'
#' @export
trimcibt <- function(x, tr = .2, alpha = .03, nboot = 2000,
                     ci.type = c("symmetric", "asymmetric"),
                     q.type = 6,
                     hyp = 0,
                     alternative = c("two.sided", "greater", "less"),
                     small.n = TRUE) {
  ci.type <- match.arg(ci.type)
  alternative <- match.arg(alternative)
  if (length(small.n) != 1L || !is.logical(small.n) || is.na(small.n)) {
    stop("small.n must be TRUE or FALSE.")
  }
  if (!is.numeric(tr) || length(tr) != 1L || !is.finite(tr) ||
      tr < 0 || tr >= 0.5) {
    stop("tr must be a single finite number in [0, 0.5).")
  }

  estimate <- mean(x, trim = tr)
  se <- trimse(x, tr = tr)
  test.stat <- (estimate - hyp) / se

  data <- matrix(
    sample(x, size = length(x) * nboot, replace = TRUE),
    nrow = nboot
  )
  data <- data - estimate
  boot.estimates <- apply(data, 1, mean, trim = tr)
  boot.t <- boot.estimates / apply(data, 1, trimse, tr = tr)

  if (alternative == "two.sided") {
    if (ci.type == "symmetric") {
      critical <- quantile(abs(boot.t), 1 - alpha, names = FALSE, type = q.type)
      ci <- c(estimate - critical * se, estimate + critical * se)
      tail.count <- sum(abs(boot.t) >= abs(test.stat))
      p.value <- if (small.n) {
        (tail.count + 1) / (nboot + 1)
      } else {
        tail.count / nboot
      }
    } else { # "asymmetric"
      lower.critical <- quantile(boot.t, alpha / 2, names = FALSE, type = q.type)
      upper.critical <- quantile(boot.t, 1 - alpha / 2, names = FALSE, type = q.type)
      ci <- c(
        estimate - upper.critical * se,
        estimate - lower.critical * se
      )
      tail.count <- if (test.stat < 0) {
        sum(test.stat > boot.t)
      } else {
        sum(test.stat < boot.t)
      }
      tail.p <- if (small.n) {
        (tail.count + 1) / (nboot + 1)
      } else {
        tail.count / nboot
      }
      p.value <- min(1, 2 * tail.p)
    }
  } else if (alternative == "greater") {
    critical <- quantile(boot.t, 1 - alpha, names = FALSE, type = q.type)
    ci <- c(estimate - critical * se, Inf)
    tail.count <- sum(boot.t >= test.stat)
    p.value <- if (small.n) {
      (tail.count + 1) / (nboot + 1)
    } else {
      tail.count / nboot
    }
  } else { # "less"
    critical <- quantile(boot.t, alpha, names = FALSE, type = q.type)
    ci <- c(-Inf, estimate - critical * se)
    tail.count <- sum(boot.t <= test.stat)
    p.value <- if (small.n) {
      (tail.count + 1) / (nboot + 1)
    } else {
      tail.count / nboot
    }
  }

  list(
    estimate = estimate,
    ci = ci,
    test.stat = test.stat,
    p.value = p.value,
    n = length(x),
    ci.type = ci.type,
    q.type = q.type,
    hyp = hyp,
    alternative = alternative,
    small.n = small.n,
    boot.estimates = boot.estimates,
    boot.t = boot.t
  )
}

#' Parametric inference for two independent trimmed means
#'
#' Compares the trimmed means of two independent samples using Yuen's method.
#'
#' @param x Numeric data vector for the first group. Missing values are not removed.
#'   A two-column matrix, data frame, or list can also be supplied when `y` is
#'   omitted.
#' @param y Optional numeric data vector for the second group. Missing values
#'   are not removed.
#' @param tr Trimming proportion in [0, 0.5), with 0.2 default.
#' @param alpha Significance level, with 0.05 default.
#' @param alternative Alternative hypothesis: "two.sided", "greater", or
#'   "less".
#'
#' @return A list containing the sample sizes, trimmed mean estimates,
#'   difference in trimmed means, confidence interval, p-value, standard
#'   error, test statistic, degrees of freedom, and alternative hypothesis.
#'
#' @details
#' When `y` is omitted, `x` can be a two-column matrix or data frame, or a
#' list containing the two groups.
#'
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
#'
#' @examples
#' set.seed(666)
#' x <- rnorm(20, mean = 1)
#' y <- rnorm(20)
#' yuen(x, y)
#' yuen(x, y, alternative = "greater")
#'
#' @seealso [yuenbt()], [yuend()]
#' @export
yuen <- function(x, y = NULL, tr = .2, alpha = .05,
                 alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  if (is.null(y)) {
    if (is.matrix(x) || is.data.frame(x)) {
      y <- x[, 2]
      x <- x[, 1]
    }
    if (is.list(x)) {
      y <- x[[2]]
      x <- x[[1]]
    }
  }
  if (!is.numeric(tr) || length(tr) != 1L || !is.finite(tr) ||
    tr < 0 || tr >= 0.5) {
    stop("tr must be a single finite number in [0, 0.5).")
  }
  h1 <- length(x) - 2 * floor(tr * length(x))
  h2 <- length(y) - 2 * floor(tr * length(y))
  q1 <- (length(x) - 1) * winvar(x, tr) / (h1 * (h1 - 1))
  q2 <- (length(y) - 1) * winvar(y, tr) / (h2 * (h2 - 1))
  df <- (q1 + q2)^2 / ((q1^2 / (h1 - 1)) + (q2^2 / (h2 - 1)))
  estimate.x <- mean(x, trim = tr)
  estimate.y <- mean(y, trim = tr)
  dif <- estimate.x - estimate.y
  se <- sqrt(q1 + q2)
  test.stat <- dif / se

  if (alternative == "two.sided") {
    crit <- qt(1 - alpha / 2, df)
    ci <- dif + c(-1, 1) * crit * se
    p.value <- 2 * pt(-abs(test.stat), df)
  } else if (alternative == "greater") {
    crit <- qt(1 - alpha, df)
    ci <- c(dif - crit * se, Inf)
    p.value <- pt(test.stat, df, lower.tail = FALSE)
  } else { # "less"
    crit <- qt(1 - alpha, df)
    ci <- c(-Inf, dif + crit * se)
    p.value <- pt(test.stat, df)
  }

  list(
    n1 = length(x),
    n2 = length(y),
    est.1 = estimate.x,
    est.2 = estimate.y,
    ci = ci,
    p.value = p.value,
    dif = dif,
    se = se,
    teststat = test.stat,
    crit = crit,
    df = df,
    alternative = alternative
  )
}

#' Bootstrap-t inference for two independent trimmed means
#'
#' @param x Numeric vector for group 1. Missing values must be removed first.
#' @param y Numeric vector for group 2. Missing values must be removed first.
#' @param tr Trimming proportion in [0, 0.5).
#' @param alpha Significance level.
#' @param nboot Number of bootstrap samples.
#' @param ci.type Type of confidence interval for two-sided inference: "symmetric" or
#'   "asymmetric". Ignored for one-sided inference, which always uses signed
#'   bootstrap statistics.
#' @param q.type Type of quantile estimator to compute bootstrap interval. Default to 6.
#' @param alternative Alternative hypothesis: "two.sided", "greater", or "less".
#' @param hyp Hypothesis value for the difference in trimmed means (default to null).
#' @param small.n Logical, default is TRUE to apply the small-sample correction
#'   `(tail.count + 1) / (nboot + 1)` to bootstrap p-values.
#'
#' @return A list containing the confidence interval, test statistic,
#'   p-value, estimates, sample sizes, bootstrap estimates, and bootstrap
#'   t-statistics.
#'
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
#'
#' @examples
#' set.seed(666)
#' x <- rnorm(20, mean = 1)
#' y <- rnorm(20, mean = 0)
#'
#' yuenbt(x, y, ci.type = "symmetric", alternative = "two.sided")
#' yuenbt(x, y, ci.type = "asymmetric", alternative = "greater")
#' yuenbt(x, y, ci.type = "asymmetric", alternative = "less")
#'
#' @export
yuenbt <- function(x, y, tr = .2, alpha = .03, nboot = 2000,
                   ci.type = c("symmetric", "asymmetric"),
                   q.type = 6,
                   alternative = c("two.sided", "greater", "less"),
                   hyp = 0, small.n = TRUE) {
  ci.type <- match.arg(ci.type)
  alternative <- match.arg(alternative)
  if (length(small.n) != 1L || !is.logical(small.n) || is.na(small.n)) {
    stop("small.n must be TRUE or FALSE.")
  }

  if (!is.numeric(tr) || length(tr) != 1L || !is.finite(tr) ||
    tr < 0 || tr >= 0.5) {
    stop("tr must be a single finite number in [0, 0.5).")
  }

  estimate.x <- mean(x, trim = tr)
  estimate.y <- mean(y, trim = tr)
  estimate.dif <- estimate.x - estimate.y

  se <- sqrt(trimse(x, tr = tr)^2 + trimse(y, tr = tr)^2)
  test.stat <- (estimate.dif - hyp) / se

  xcen <- x - estimate.x # trimmmed mean centre the samples
  ycen <- y - estimate.y

  datax <- matrix( # bootstrap samples
    sample(xcen, length(x) * nboot, replace = TRUE),
    nrow = nboot
  )
  datay <- matrix(
    sample(ycen, length(y) * nboot, replace = TRUE),
    nrow = nboot
  )

  boot.estimates <- apply(datax, 1, mean, trim = tr) -
    apply(datay, 1, mean, trim = tr)
  tboot <- boot.estimates / sqrt( # bootstrap t-values
    apply(datax, 1, trimse, tr = tr)^2 +
      apply(datay, 1, trimse, tr = tr)^2
  )

  if (alternative == "two.sided") {
    if (ci.type == "symmetric") {
      critical <- quantile(abs(tboot), 1 - alpha, names = FALSE, type = q.type)
      ci <- estimate.dif + c(-1, 1) * critical * se
      tail.count <- sum(abs(tboot) >= abs(test.stat))
      p.value <- if (small.n) {
        (tail.count + 1) / (nboot + 1)
      } else {
        tail.count / nboot
      }
    } else { # "asymmetric"
      lower.critical <- quantile(tboot, alpha / 2, names = FALSE, type = q.type)
      upper.critical <- quantile(tboot, 1 - alpha / 2, names = FALSE, type = q.type)
      ci <- c(
        estimate.dif - upper.critical * se,
        estimate.dif - lower.critical * se
      )
      tail.count <- if (test.stat < 0) {
        sum(test.stat > tboot)
      } else {
        sum(test.stat < tboot)
      }
      tail.p <- if (small.n) {
        (tail.count + 1) / (nboot + 1)
      } else {
        tail.count / nboot
      }
      p.value <- min(1, 2 * tail.p)
    }
  } else if (alternative == "greater") {
    critical <- quantile(tboot, 1 - alpha, names = FALSE, type = q.type)
    ci <- c(estimate.dif - critical * se, Inf)
    tail.count <- sum(tboot >= test.stat)
    p.value <- if (small.n) {
      (tail.count + 1) / (nboot + 1)
    } else {
      tail.count / nboot
    }
  } else { # "less"
    critical <- quantile(tboot, alpha, names = FALSE, type = q.type)
    ci <- c(-Inf, estimate.dif - critical * se)
    tail.count <- sum(tboot <= test.stat)
    p.value <- if (small.n) {
      (tail.count + 1) / (nboot + 1)
    } else {
      tail.count / nboot
    }
  }

  list(
    ci = ci,
    test.stat = test.stat,
    p.value = p.value,
    est.1 = estimate.x,
    est.2 = estimate.y,
    estimate = estimate.dif,
    se = se,
    ci.type = ci.type,
    q.type = q.type,
    alternative = alternative,
    small.n = small.n,
    boot.estimates = boot.estimates,
    boot.t = tboot,
    n1 = length(x),
    n2 = length(y)
  )
}

#' Parametric inference for two dependent trimmed means
#'
#' Compares the trimmed means of two paired samples using Yuen's method.
#'
#' @param x Numeric data vector for the first measurement.
#' @param y Numeric data vector for the second measurement. Must have the same
#'   length as `x`.
#' @param tr Trimming proportion in [0, 0.5), with 0.2 default.
#' @param alpha Significance level, with 0.05 default.
#'
#' @return A list containing the confidence interval, p-value, trimmed mean
#'   estimates, difference in trimmed means, standard error, test statistic,
#'   sample size, and degrees of freedom.
#'
#' @details
#' Missing values must be removed before calling this function. For inferences based on difference scores, use trimci() or trimcibt().
#'
#' @export
#'
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
#'
#' @examples
#' set.seed(666)
#' x <- rnorm(20, mean = 1)
#' y <- x + rnorm(20)
#' yuend(x, y)
#'
#' @seealso [trimci()], [trimcibt()], [yuen()]
yuend <- function(x, y, tr = .2, alpha = .05) {
  if (length(x) != length(y)) stop("x and y must have equal sample sizes.")
  h1 <- length(x) - 2 * floor(tr * length(x))
  q1 <- (length(x) - 1) * winvar(x, tr)
  q2 <- (length(y) - 1) * winvar(y, tr)
  q3 <- (length(x) - 1) * wincor(x, y, tr)$cov
  df <- h1 - 1
  se <- sqrt((q1 + q2 - 2 * q3) / (h1 * (h1 - 1)))
  crit <- qt(1 - alpha / 2, df)
  dif <- mean(x, tr) - mean(y, tr)
  low <- dif - crit * se
  up <- dif + crit * se
  test <- dif / se
  p.value <- 2 * (1 - pt(abs(test), df))
  list(
    ci = c(low, up),
    p.value = p.value,
    est1 = mean(x, tr),
    est2 = mean(y, tr),
    dif = dif,
    se = se,
    teststat = test,
    n = length(x),
    df = df
  )
}

#' Bootstrap-t inference for two dependent trimmed means
#'
#' Computes a bootstrap-t confidence interval and p-value for the difference between the
#' marginal trimmed means of paired data.
#'
#' @param x Numeric data vector for the first measurement.
#' @param y Numeric data vector for the second measurement. Must have the same
#'   length as `x`.
#' @param tr Trimming proportion in [0, 0.5), with 0.2 default.
#' @param alpha Significance level, with 0.05 default.
#' @param nboot Number of bootstrap samples, with 2000 default.
#' @param ci.type Type of confidence interval for two-sided inference: "symmetric" or
#'   "asymmetric". Ignored for one-sided inference, which always uses signed
#'   bootstrap statistics.
#' @param q.type Type of quantile estimator to compute the bootstrap interval.
#'   Default to 6.
#' @param alternative Alternative hypothesis: "two.sided", "greater", or "less".
#' @param hyp Hypothesis value for the difference in trimmed means. Default to 0.
#' @param small.n Default to TRUE to apply the small-sample correction
#'   `(tail.count + 1) / (nboot + 1)` to the bootstrap p-value.
#'
#' @return A list containing the confidence interval, trimmed mean estimates,
#'   difference in trimmed means, and p-value.
#'
#' @details
#' Missing values must be removed before calling this function. When making inferences on means,
#' the difference between means is the same as the mean of the difference scores. However, this
#' is not the case for trimmed means, medians and other quantiles. For inferences on the trimmed mean of difference scores,
#' use `trimcibt()` or `trimci()` for a non-bootstrap approach.
#'
#' @export
#'
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
#'
#' @examples
#' set.seed(666)
#' x <- rnorm(20, mean = 1)
#' y <- x + rnorm(20)
#' yuendbt(x, y, ci.type = "symmetric", alternative = "two.sided")
#' yuendbt(x, y, ci.type = "asymmetric", alternative = "two.sided")
#' yuendbt(x, y, alternative = "greater")
#' yuendbt(x, y, alternative = "less")
yuendbt <- function(x, y, tr = .2, alpha = .05, nboot = 2000,
                    ci.type = c("symmetric", "asymmetric"),
                    q.type = 6,
                    alternative = c("two.sided", "greater", "less"),
                    hyp = 0, small.n = TRUE) {
  # Original name was ydbt. Renamed for consistency with naming convention
  # for inferences on independent group trimmed means.
  ci.type <- match.arg(ci.type)
  alternative <- match.arg(alternative)
  if (length(small.n) != 1L || !is.logical(small.n) || is.na(small.n)) {
    stop("small.n must be TRUE or FALSE.")
  }
  if (length(x) != length(y)) stop("x and y must have equal sample sizes.")
  data <- matrix(sample(length(y), size = length(y) * nboot, replace = TRUE), nrow = nboot)
  e1 <- mean(x, tr)
  e2 <- mean(y, tr)
  xcen <- x - mean(x, tr)
  ycen <- y - mean(y, tr)
  bvec <- apply(data, 1, tsub, xcen, ycen, tr)
  # bvec is a 1 by nboot matrix containing the bootstrap test statistics.
  dotest <- yuend(x, y, tr = tr)
  estse <- dotest$se
  dif <- mean(x, tr) - mean(y, tr)
  test.stat <- (dif - hyp) / estse

  if (alternative == "two.sided") {
    if (ci.type == "symmetric") {
      critical <- quantile(abs(bvec), 1 - alpha, names = FALSE, type = q.type)
      ci <- dif + c(-1, 1) * critical * estse
      tail.count <- sum(abs(bvec) >= abs(test.stat))
      p.value <- if (small.n) {
        (tail.count + 1) / (nboot + 1)
      } else {
        tail.count / nboot
      }
    } else {
      lower.critical <- quantile(bvec, alpha / 2, names = FALSE, type = q.type)
      upper.critical <- quantile(bvec, 1 - alpha / 2, names = FALSE, type = q.type)
      ci <- c(dif - upper.critical * estse, dif - lower.critical * estse)
      tail.count <- if (test.stat < 0) {
        sum(test.stat > bvec)
      } else {
        sum(test.stat < bvec)
      }
      tail.p <- if (small.n) {
        (tail.count + 1) / (nboot + 1)
      } else {
        tail.count / nboot
      }
      p.value <- min(1, 2 * tail.p)
    }
  } else if (alternative == "greater") {
    critical <- quantile(bvec, 1 - alpha, names = FALSE, type = q.type)
    ci <- c(dif - critical * estse, Inf)
    tail.count <- sum(bvec >= test.stat)
    p.value <- if (small.n) {
      (tail.count + 1) / (nboot + 1)
    } else {
      tail.count / nboot
    }
  } else {
    critical <- quantile(bvec, alpha, names = FALSE, type = q.type)
    ci <- c(-Inf, dif - critical * estse)
    tail.count <- sum(bvec <= test.stat)
    p.value <- if (small.n) {
      (tail.count + 1) / (nboot + 1)
    } else {
      tail.count / nboot
    }
  }
  #bsort <- if (ci.type == "symmetric") sort(abs(bvec)) else sort(bvec)
  list(
    ci = ci,
    Est.1 = e1,
    Est.2 = e2,
    dif = dif,
    p.value = p.value,
    test.stat = test.stat,
    ci.type = ci.type,
    q.type = q.type,
    alternative = alternative,
    hyp = hyp,
    small.n = small.n,
    boot.t = bvec
  )
}

#' Compute a bootstrap test statistic for paired trimmed means
#'
#' @param isub Bootstrap sample indices.
#' @param x Numeric data vector for the first measurement.
#' @param y Numeric data vector for the second measurement.
#' @param tr Trimming proportion.
#' @return A bootstrap test statistic.
#' @keywords internal
tsub <- function(isub,x,y,tr){
  tsub <- yuend(x[isub],y[isub],tr=tr)$teststat
  tsub
}

#' Winsorized correlation and covariance
#'
#' Computes the Winsorized correlation and covariance between variables.
#' Missing values must be removed before calling this function.
#'
#' @param x Numeric vector, or a numeric matrix with at least two columns.
#' @param y Optional numeric vector. If supplied, it is paired with `x`.
#' @param tr Winsorization proportion in [0, 0.5), with 0.2 default.
#' @return A list containing the sample size, Winsorized correlation,
#'   Winsorized covariance, and p-value.
#'
#' @export
wincor <- function(x, y = NULL, tr = .2) {
  data <- if (is.null(y)) x else cbind(x, y)
  if (!is.matrix(data)) {
    stop("The data must be stored in a matrix when y is NULL.")
  }
  nval <- nrow(data)

  if (ncol(data) == 2) {
    result <- wincor.sub(data[, 1], data[, 2], tr)
    return(list(
      n = nval,
      cor = result$cor,
      cov = result$cov,
      p.value = result$p.value
    ))
  }

  if (ncol(data) < 2) {
    stop("The data must contain at least two columns.")
  }
  wcor <- matrix(1, ncol(data), ncol(data))
  wcov <- matrix(0, ncol(data), ncol(data))
  siglevel <- matrix(NA_real_, ncol(data), ncol(data))
  for (i in seq_len(ncol(data))) {
    for (j in i:ncol(data)) {
      result <- wincor.sub(data[, i], data[, j], tr)
      wcor[i, j] <- result$cor
      wcor[j, i] <- wcor[i, j]
      wcov[i, j] <- result$cov
      wcov[j, i] <- wcov[i, j]
      if (i != j) {
        siglevel[i, j] <- result$p.value
        siglevel[j, i] <- siglevel[i, j]
      }
    }
  }
  list(n = nval, cor = wcor, cov = wcov, p.value = siglevel)
}

#' Compute a Winsorized correlation for two numeric vectors
#'
#' @param x Numeric vector.
#' @param y Numeric vector of the same length as `x`.
#' @param tr Winsorization proportion in [0, 0.5), with 0.2 default.
#' @return A list containing the Winsorized correlation, covariance, and p-value.
#' @keywords internal
wincor.sub <- function(x, y, tr = .2) {
  n <- length(x)
  g <- floor(tr * n)
  xsort <- sort(x)
  ysort <- sort(y)
  xvec <- pmin(pmax(x, xsort[g + 1]), xsort[n - g])
  yvec <- pmin(pmax(y, ysort[g + 1]), ysort[n - g])
  wcor <- cor(xvec, yvec)
  wcov <- var(xvec, yvec)
  p.value <- NA_real_
  if (!all(x == y)) {
    test <- wcor * sqrt((n - 2) / (1 - wcor^2))
    p.value <- 2 * (1 - pt(abs(test), n - 2 * g - 2))
  }
  list(cor = wcor, cov = wcov, p.value = p.value)
}

#' Parametric standard error of a trimmed mean
#'
#' @param x Numeric data vector. Missing values must be removed first.
#' @param tr Trimming proportion in [0, 0.5), with default to 0.2.
#' @return The estimated standard error of the trimmed mean.
#' @export
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
trimse <- function(x, tr = .2) {
  if (!is.numeric(tr) || length(tr) != 1L || !is.finite(tr) || tr < 0 || tr >= 0.5) {
    stop("tr must be a single finite number in [0, 0.5).")
  }
  n <- length(x)
  trimse <- sqrt(winvar(x, tr)) / ((1 - 2 * tr) * sqrt(n))
  trimse
}

#' Winsorized variance
#'
#' @param x Numeric data vector. Missing values must be removed first.
#' @param tr Winsorization proportion in [0, 0.5).
#' @return The Winsorized sample variance.
#' @export
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
winvar <- function(x, tr = .2) {
  if (!is.numeric(tr) || length(tr) != 1L || !is.finite(tr) || tr < 0 || tr >= 0.5) {
    stop("tr must be a single finite number in [0, 0.5).")
  }
  remx <- x
  if (any(is.na(remx))) {
    wv <- NA
  } else {
    y <- sort(x)
    n <- length(x)
    ibot <- floor(tr * n) + 1
    itop <- n - ibot + 1
    xbot <- y[ibot]
    xtop <- y[itop]
    y <- ifelse(y <= xbot, xbot, y)
    y <- ifelse(y >= xtop, xtop, y)
    wv <- var(y)
  }
  wv
}

#' Parametric confidence interval for a median
#'
#' Uses the Hettmansperger-Sheather interpolation method.
#'
#' @param x Numeric data vector. Missing values must be removed first.
#' @param alpha Significance level, with default to 0.03 (97% confidence).
#' @return A length-two numeric confidence interval.
#'
#' @details
#' If there are duplicate values in `x`, a gain in power might be obtained by calling `onesampb(x, est=hd)` instead.
#' This will provide a percentile bootstrap confidence interval for the Harrell-Davis estimator of the median.
#'
#' @export
#'
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
sint <- function(x, alpha = .03) {
  if (any(duplicated(x))) {
    warning("Duplicate values detected: onesampb(x, est=hd) might have more power")
  }
  n <- length(x)
  k <- qbinom(alpha / 2, n, .5)
  gk <- pbinom(n - k, n, .5) - pbinom(k - 1, n, .5)
  if (gk >= 1 - alpha) {
    gkp1 <- pbinom(n - k - 1, n, .5) - pbinom(k, n, .5)
    kp <- k + 1
  } else {
    k <- k - 1
    gk <- pbinom(n - k, n, .5) - pbinom(k - 1, n, .5)
    gkp1 <- pbinom(n - k - 1, n, .5) - pbinom(k, n, .5)
    kp <- k + 1
  }
  xsort <- sort(x)
  nmk <- n - k
  nmkp <- nmk + 1
  ival <- (gk - 1 + alpha) / (gk - gkp1)
  lam <- ((n - k) * ival) / (k + (n - 2 * k) * ival)
  low <- lam * xsort[kp] + (1 - lam) * xsort[k]
  hi <- lam * xsort[nmk] + (1 - lam) * xsort[nmkp]
  ci <- c(low, hi)
  ci
}

#' Percentile bootstrap confidence interval and p-value
#'
#' Standalone utility to (re)compute a percentile bootstrap confidence
#' interval and p-value from bootstrap estimates, for any combination of
#' `hyp`, `alternative`, and `alpha`. This function is used
#' internally by [onesampb()] and [twosampb()], so it can be applied to the
#' result list returned by either of those functions without re-running the
#' bootstrap.
#'
#' @param x A numeric vector of bootstrap estimates, or a list of results
#'   from [onesampb()] or [twosampb()] (or any list containing a numeric
#'   `boot.estimates` element).
#' @param alpha Significance level, with default to 0.03 (97% intended coverage).
#' @param hyp Hypothesis value, with default to 0.
#' @param alternative Alternative hypothesis: "two.sided", "greater", or "less".
#' @param q.type Type of quantile estimator. Default to 6.
#' @param small.n Logical; apply the small-sample correction
#'   `(count + 1) / (nboot + 1)` to the bootstrap p-value.
#' @details
#' For "greater" the confidence interval is `[lower bound, Inf]`.
#' For "less" the confidence interval is `[-Inf, upper bound]`.
#' @return A list containing the confidence interval, p-value, hyp,
#'   alternative, and alpha.
#' @export
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
pbci <- function(x, alpha = .03, hyp = 0,
                alternative = c("two.sided", "greater", "less"),
                q.type = 6, small.n = TRUE) {
  alternative <- match.arg(alternative)
  if (length(small.n) != 1L || !is.logical(small.n) || is.na(small.n)) {
    stop("small.n must be TRUE or FALSE.")
  }
  boot.estimates <- if (is.list(x)) x$boot.estimates else x
  if (!is.numeric(boot.estimates)) {
    stop("x must be a numeric vector, or a list containing a numeric boot.estimates element.")
  }
  nboot <- length(boot.estimates)
  tie.count <- sum(boot.estimates == hyp)
  upper.count <- sum(boot.estimates > hyp) + 0.5 * tie.count
  lower.count <- sum(boot.estimates < hyp) + 0.5 * tie.count
  upper.p <- if (small.n) {
    (upper.count + 1) / (nboot + 1)
  } else {
    upper.count / nboot
  }
  lower.p <- if (small.n) {
    (lower.count + 1) / (nboot + 1)
  } else {
    lower.count / nboot
  }
  if (alternative == "two.sided") {
    ci <- quantile(boot.estimates, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE, type = q.type)
    p.value <- min(1, 2 * min(upper.p, lower.p))
  } else if (alternative == "greater") {
    ci <- c(quantile(boot.estimates, alpha, names = FALSE, type = q.type), Inf)
    p.value <- lower.p
  } else {
    ci <- c(-Inf, quantile(boot.estimates, 1 - alpha, names = FALSE, type = q.type))
    p.value <- upper.p
  }
  list(ci = ci, p.value = p.value, hyp = hyp, alternative = alternative, alpha = alpha, small.n = small.n)
}

#' One-sample percentile bootstrap inference
#'
#' @param x Numeric data vector. Missing values must be removed first.
#' @param est Estimator function, with default to hd().
#' @param alpha Significance level, with default to 0.03 (97% intended coverage).
#' @param nboot Number of bootstrap samples, with default to 2000.
#' @param hyp Hypothesis value, with default to 0.
#' @param alternative Alternative hypothesis: "two.sided", "greater", or "less".
#' @param q.type Type of quantile estimator. Default to 6.
#' @param small.n Logical; apply the small-sample correction
#'   `(count + 1) / (nboot + 1)` to the bootstrap p-value.
#' @param ... Additional arguments passed to `est`, for instance `tr` for the trimmed mean or `q` for `hd`.
#' @return A list containing the percentile confidence interval, estimate,
#'   p-value, sample size, and bootstrap estimates.
#' @export
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
onesampb <- function(x, est = hd, alpha = .03, nboot = 2000, hyp = 0,
                     alternative = c("two.sided", "greater", "less"),
                     q.type = 6, small.n = TRUE, ...) {
  alternative <- match.arg(alternative)
  if (length(small.n) != 1L || !is.logical(small.n) || is.na(small.n)) {
    stop("small.n must be TRUE or FALSE.")
  }
  n <- length(x)
  estimate <- est(x, ...)
  data <- matrix(sample(x, size = n * nboot, replace = TRUE), nrow = nboot)
  boot.estimates <- apply(data, 1, est, ...)
  result <- pbci(boot.estimates, alpha, hyp, alternative, q.type, small.n)

  list(
    ci = result$ci,
    n = n,
    estimate = estimate,
    p.value = result$p.value,
    hyp = hyp,
    alternative = alternative,
    small.n = small.n,
    boot.estimates = boot.estimates
  )
}

#' Two-sample percentile bootstrap inference
#'
#' Compares an estimator between independent groups or paired observations.
#' Missing values are not removed internally.
#'
#' @param x Numeric vector for group 1, or a two-column numeric matrix when
#'   `ind = FALSE` and `y` is omitted.
#' @param y Numeric vector for group 2. Omit when `x` is a two-column matrix.
#' @param alpha Significance level, with default to 0.03 (97% intended coverage).
#' @param nboot Number of bootstrap samples.
#' @param est Estimator function. Defaults to [hd()].
#' @param hyp Null hypothesis value for the difference.
#' @param alternative Alternative hypothesis: "two.sided", "greater", or "less".
#' @param q.type Type of quantile estimator. Default to 6.
#' @param ind Logical; if `TRUE` (default), resample the two groups
#'   independently. If `FALSE`, resample rows of paired observations; the two groups must have
#'   equal lengths.
#' @param small.n Logical; apply the small-sample correction
#'   `(count + 1) / (nboot + 1)` to the bootstrap p-value.
#' @param ... Additional arguments passed to `est`.
#' @return A list containing estimates, confidence interval, p-value, sample
#'   sizes, bootstrap estimate differences, and bootstrap variance.
#' @export
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
twosampb <- function(x, y = NULL, alpha = .03, nboot = 2000, est = hd,
                     hyp = 0, alternative = c("two.sided", "greater", "less"),
                     q.type = 6, ind = TRUE, small.n = TRUE, ...) {
  alternative <- match.arg(alternative)
  if (length(small.n) != 1L || !is.logical(small.n) || is.na(small.n)) {
    stop("small.n must be TRUE or FALSE.")
  }
  if (!ind && is.null(y)) {
    if (!is.matrix(x) || ncol(x) != 2L) {
      stop("For ind = FALSE with y omitted, x must be a matrix with two columns.")
    }
    y <- x[, 2]
    x <- x[, 1]
  }
  if (is.null(y)) {
    stop("y must be supplied unless x is a two-column matrix and ind = FALSE.")
  }
  if (!ind && length(x) != length(y)) {
    stop("Paired x and y must have the same length.")
  }

  if (!ind) {
    indices <- matrix(
      sample(seq_along(x), length(x) * nboot, replace = TRUE),
      nrow = nboot
    )
    boot.estimates.x <- apply(indices, 1, function(index) {
      est(x[index], ...)
    })
    boot.estimates.y <- apply(indices, 1, function(index) {
      est(y[index], ...)
    })
    boot.estimates <- boot.estimates.x - boot.estimates.y
  } else {
    datax <- matrix(
      sample(x, length(x) * nboot, replace = TRUE),
      nrow = nboot
    )
    datay <- matrix(
      sample(y, length(y) * nboot, replace = TRUE),
      nrow = nboot
    )
    boot.estimates.x <- apply(datax, 1, est, ...)
    boot.estimates.y <- apply(datay, 1, est, ...)
    boot.estimates <- boot.estimates.x - boot.estimates.y
  }

  result <- pbci(boot.estimates, alpha, hyp, alternative, q.type, small.n)
  estimate.x <- est(x, ...)
  estimate.y <- est(y, ...)
  list(
    est.1 = estimate.x,
    est.2 = estimate.y,
    estimate = estimate.x - estimate.y,
    ci = result$ci,
    p.value = result$p.value,
    hyp = hyp,
    alternative = alternative,
    small.n = small.n,
    sq.se = var(boot.estimates),
    boot.estimates = boot.estimates,
    boot.estimates.x = boot.estimates.x,
    boot.estimates.y = boot.estimates.y,
    ind = ind,
    n1 = length(x),
    n2 = length(y)
  )
}

#' Evaluate Huber's Psi function for each value in the vector x
#'
#' @param x Numeric values.
#' @param bend Positive bending constant, with default to 1.28.
#' @return The Huber Psi transformation of `x`.
#' @export
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
hpsi <- function(x, bend = 1.28) {
  hpsi <- ifelse(abs(x) <= bend, x, bend * sign(x))
  hpsi
}

#' One-step M-estimator of location
#'
#' @param x Numeric data vector.
#' @param bend Huber's Psi bending constant, with default to 1.28.
#' @param MED Logical, with default to TRUE to use the median as the initial location estimate, otherwise use modified one-step M-estimator.
#' @return A one-step M-estimate of location.
#' @export
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
onestep <- function(x, bend = 1.28, MED = TRUE) {
  #
  #  Compute one-step M-estimator of location using Huber's Psi.
  #  The default bending constant is 1.28
  #
  #  MED=TRUE: initial estimate is the median
  #  Otherwise use modified one-step M-estimator
  #
  if (MED) init.loc <- median(x)
  if (!MED) init.loc <- mom(x, bend = bend)
  y <- (x - init.loc) / mad(x) # mad in splus is madn in the book.
  A <- sum(hpsi(y, bend))
  B <- length(x[abs(y) <= bend])
  onestep <- median(x) + mad(x) * A / B
  onestep
}

#' Modified one-step MOM location estimator
#'
#' @param x Numeric data vector.
#' @param bend Bending constant.
#' @return A MOM estimate of location.
#' @export
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
mom <- function(x, bend = 2.24) {
  flag1 <- (x > median(x) + bend * mad(x))
  flag2 <- (x < median(x) - bend * mad(x))
  flag <- rep(TRUE, length(x))
  flag[flag1] <- FALSE
  flag[flag2] <- FALSE
  mom <- mean(x[flag])
  mom
}

#' Evaluate an empirical cumulative distribution function
#'
#' Computes the proportion of values in `x` less than or equal to `val`.
#' This internal helper is used by `ks()`.
#'
#' @param x Numeric data vector.
#' @param val Value at which to evaluate the empirical distribution function.
#' @return The proportion of values in `x` less than or equal to `val`.
#' @keywords internal
ks_ecdf <- function(x, val) {
  sum(x <= val) / length(x)
}

#' Kolmogorov-Smirnov test for two samples
#'
#' Computes unweighted or weighted Kolmogorov-Smirnov test statistics and
#' their critical values. For tied data, the reported p-value is exact only
#' for the unweighted test.
#'
#' @param x Numeric data vector for the first group.
#' @param y Numeric data vector for the second group.
#' @param w Logical; compute the weighted statistic.
#' @param sig Logical; compute the significance level.
#' @param alpha Significance level for the critical value.
#' @return A list containing the test statistic, critical value, and p-value.
#' @export
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
ks <- function(x,y,w=FALSE,sig=TRUE,alpha=.05){
  y<-y[!is.na(y)]
n1 <- length(x)
n2 <- length(y)
w<-as.logical(w)
sig<-as.logical(sig)
tie<-logical(1)
siglevel<-NA
z<-sort(c(x,y))  # Pool and sort the observations
tie=FALSE
chk=sum(duplicated(x,y))
if(chk>0)tie=TRUE
v<-1   # Initializes v
for (i in 1:length(z))v[i]<-abs(ks_ecdf(x,z[i]) - ks_ecdf(y,z[i]))
ks<-max(v)
if(!tie)crit=ks.crit(n1=n1,n2=n2,alpha=alpha)
else crit=ksties.crit(x,y,alpha=alpha)
if(!w && sig && !tie)siglevel<-kssig(length(x),length(y),ks)
if(!w && sig && tie)siglevel<-kstiesig(x,y,ks)
if(w){
crit=ksw.crit(length(x),length(y),alpha=alpha)
for (i in 1:length(z)){
temp<-(length(x)*ks_ecdf(x,z[i])+length(y)*ks_ecdf(y,z[i]))/length(z)
temp<-temp*(1.-temp)
v[i]<-v[i]/sqrt(temp)
}
v<-v[!is.na(v)]
ks<-max(v)*sqrt(length(x)*length(y)/length(z))
if(sig)siglevel<-kswsig(length(x),length(y),ks)
if(tie && sig)
warning(paste("Ties were detected. The reported significance level of the
weighted Kolmogorov-Smirnov test statistic is not exact."))
}
list(test=ks,critval=crit,p.value=siglevel)
}

# -----------------------------------------------------------------------------
# ks subfunctions -- do not export
ks.crit<-function(n1,n2,alpha=.05){
  #
  # Compute a critical value so that probability coverage is approximately
  # 1-alpha
  #
  START=sqrt(0-log(alpha/2)*(n1+n2)/(2*n1*n2))
  crit=optim(START,ks.sub,n1=n1,n2=n2,alpha=alpha,lower=.001,upper=.86,method='Brent')$par
  crit
}

ksties.crit<-function(x,y,alpha=.05){
#
# Compute a critical value so that probability coverage is approximately
# 1-alpha
#
n1=length(x)
n2=length(y)
START=sqrt(0-log(alpha/2)*(n1+n2)/(2*n1*n2))
crit=optim(START,ksties.sub,x=x,y=y,alpha=alpha,lower=.001,upper=.86,method='Brent')$par
crit
}

ksties.sub<-function(crit,x,y,alpha){
  v=kstiesig(x,y,crit)
  dif=abs(alpha-v)
  dif
}

ks.sub<-function(crit,n1,n2,alpha){
  v=kssig(n1,n2,crit)
  dif=abs(alpha-v)
  dif
}

ksw.crit<-function(n1,n2,alpha=.05){
  #
  # Compute a critical value so that probability coverage is
  # >= 1-alpha while being close as possible to 1-alpha
  #
  if(alpha>.1)stop('The function assumes alpha is at least .1')
  crit=2.4
  del=.05
  pc=.12
  while(pc>alpha){
    crit=crit+.05
    pc=kswsig(n1,n2,crit)
  }
  crit
}

kstiesig<-function(x,y,val){
  #
  #    Compute significance level of the  Kolmogorov-Smirnov test statistic
  #    for the data in x and y.
  #    This function allows ties among the  values.
  #    val=observed value of test statistic
  #
  m<-length(x)
  n<-length(y)
  z<-c(x,y)
  z<-sort(z)
  cmat<-matrix(0,m+1,n+1)
  umat<-matrix(0,m+1,n+1)
  for (i in 0:m){
    for (j in 0:n){
      if(abs(i/m-j/n)<=val)cmat[i+1,j+1]<-1e0
      k<-i+j
      if(k > 0 && k<length(z) && z[k]==z[k+1])cmat[i+1,j+1]<-1
    }
  }
  for (i in 0:m){
    for (j in 0:n)if(i*j==0)umat[i+1,j+1]<-cmat[i+1,j+1]
    else umat[i+1,j+1]<-cmat[i+1,j+1]*(umat[i+1,j]+umat[i,j+1])
  }
  term<-lgamma(m+n+1)-lgamma(m+1)-lgamma(n+1)
  kstiesig<-1.-umat[m+1,n+1]/exp(term)
  kstiesig
}

kssig<-function(m,n,val){
  #
  #    Compute significance level of the  Kolmogorov-Smirnov test statistic
  #    m=sample size of first group
  #    n=sample size of second group
  #    val=observed value of test statistic
  #
  cmat<-matrix(0,m+1,n+1)
  umat<-matrix(0,m+1,n+1)
  for (i in 0:m){
    for (j in 0:n)cmat[i+1,j+1]<-abs(i/m-j/n)
  }
  cmat<-ifelse(cmat<=val,1e0,0e0)
  for (i in 0:m){
    for (j in 0:n)if(i*j==0)umat[i+1,j+1]<-cmat[i+1,j+1]
    else umat[i+1,j+1]<-cmat[i+1,j+1]*(umat[i+1,j]+umat[i,j+1])
  }
  term<-lgamma(m+n+1)-lgamma(m+1)-lgamma(n+1)
  kssig<-1.-umat[m+1,n+1]/exp(term)
  kssig=max(0,kssig)
  kssig
}

kswsig<-function(m,n,val){
  #
  #    Compute significance level of the weighted
  #    Kolmogorov-Smirnov test statistic
  #
  #    m=sample size of first group
  #    n=sample size of second group
  #    val=observed value of test statistic
  #
  mpn<-m+n
  cmat<-matrix(0,m+1,n+1)
  umat<-matrix(0,m+1,n+1)
  for (i in 1:m-1){
    for (j in 1:n-1)cmat[i+1,j+1]<-abs(i/m-j/n)*sqrt(m*n/((i+j)*(1-(i+j)/mpn)))
  }
  cmat<-ifelse(cmat<=val,1,0)
  for (i in 0:m){
    for (j in 0:n)if(i*j==0)umat[i+1,j+1]<-cmat[i+1,j+1]
    else umat[i+1,j+1]<-cmat[i+1,j+1]*(umat[i+1,j]+umat[i,j+1])
  }
  term<-lgamma(m+n+1)-lgamma(m+1)-lgamma(n+1)
  kswsig<-1.-umat[m+1,n+1]/exp(term)
  kswsig
}
# -----------------------------------------------------------------------------

#' Calculate P(X<Y) for paired observations
#'
#' Estimates the probability that an observation in `x` is less than its paired
#' observation in `y`, with a confidence interval based on Pratt's method.
#'
#' @param x A numeric vector.
#' @param y A numeric vector with the same length as `x`.
#' @param alpha Significance level for the (1-alpha)confidence interval, with 0.05
#'   default.
#'
#' @return A list with the estimated probability (`phat`) and its confidence
#'   interval (`ci`).
#'
#' @details The vectors `x` and `y` contain dependent, paired observations.
#'
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
#'
#' @examples
#' set.seed(666)
#' x <- rnorm(20)
#' y <- x + rnorm(20)
#' pxlyd(x, y)
#'
#' @seealso [pxly()]
#' @export
pxlyd <- function(x,y,alpha=.05){
  #
  # For two dependent variables, x and y,
  # estimate p=P(X<Y)
  #
  dif<-(x<y)
  temp<-binomci(y=dif,alpha=alpha)
  phat<-temp$phat
  ci<-temp$ci
  list(phat=phat,ci=ci)
}

#' Pratt confidence interval for a binomial probability
#'
#' Computes a confidence interval for a binomial probability using Pratt's
#' method. This helper is used by [pxlyd()].
#'
#' @param x Number of observed successes. Defaults to the sum of `y`.
#' @param nn Number of trials. Defaults to the length of `y`.
#' @param y Optional numeric or logical vector of binary outcomes.
#' @param alpha Significance level, with 0.05 default.
#'
#' @return A list containing the estimated probability (`phat`), confidence
#'   interval (`ci`), and number of trials (`n`).
#'
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
#'
#' @keywords internal
binomci <- function(x = sum(y), nn = length(y), y = NULL, alpha = .05) {
  if (!is.null(y)) {
    nn <- length(y)
  }
  if (nn == 1) stop("Something is wrong: number of observations is only 1")
  n <- nn
  if (x != n && x != 0) {
    z <- qnorm(1 - alpha / 2)
    ratio <- ((x + 1) / (n - x))^2
    term.b <- 81 * (x + 1) * (n - x) - 9 * n - 8
    term.c <- -3 * z * sqrt(9 * (x + 1) * (n - x) * (9 * n + 5 - z^2) + n + 1)
    term.d <- 81 * (x + 1)^2 - 9 * (x + 1) * (2 + z^2) + 1
    denominator <- 1 + ratio * ((term.b + term.c) / term.d)^3
    upper <- 1 / denominator
    ratio <- (x / (n - x - 1))^2
    term.b <- 81 * x * (n - x - 1) - 9 * n - 8
    term.c <- 3 * z * sqrt(9 * x * (n - x - 1) * (9 * n + 5 - z^2) + n + 1)
    term.d <- 81 * x^2 - 9 * x * (2 + z^2) + 1
    denominator <- 1 + ratio * ((term.b + term.c) / term.d)^3
    lower <- 1 / denominator
  }
  if (x == 0) {
    lower <- 0
    upper <- 1 - alpha^(1 / n)
  }
  if (x == 1) {
    upper <- 1 - (alpha / 2)^(1 / n)
    lower <- 1 - (1 - alpha / 2)^(1 / n)
  }
  if (x == n - 1) {
    lower <- (alpha / 2)^(1 / n)
    upper <- (1 - alpha / 2)^(1 / n)
  }
  if (x == n) {
    lower <- alpha^(1 / n)
    upper <- 1
  }
  phat <- x / n
  list(phat = phat, ci = c(lower, upper), n = n)
}

#' Calculate P(X<Y), the probability that a random observation from vector x is less than a random observation from vector y.
#' @param x A numeric vector
#' @param y A numeric vector
#' @return A probability `[0,1]`.
#' @details The two vectors x and y are independent.
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
#'
#' Birnbaum, Z. W. (1956). On a Use of the Mann-Whitney Statistic.
#' In Proceedings of the Third Berkeley Symposium on Mathematical Statistics and Probability, Volume 1:
#' Contributions to the Theory of Statistics: 3.1 (pp. 13–18).
#' University of California Press.
#' @seealso [mxmy()]
#' @export
pxly <- function(x,y){
  est <- mean(outer(x, y, FUN="<"))
  est
}

#' Estimate the median of the distribution of all pairwise differences.
#' @param x A numeric vector
#' @param y A numeric vector
#' @param est An estimator applied to the distribution X-Y, with default to the median.
#' @param ... Additional arguments passed to `est`.
#' @return A length-one object.
#' @details The two vectors x and y are independent. The median is returned by default, but another estimator can be used.
#' @references Wilcox, Rand R. 2022. Introduction to Robust Estimation and Hypothesis Testing. 5th edn.
#' Statistical Modeling and Decision Science. San Diego, CA: Academic Press.
#'
#' Hodges, J. L., & Lehmann, E. L. (1963). Estimates of Location Based on Rank Tests.
#' The Annals of Mathematical Statistics, 34(2), 598–611.
#' https://doi.org/10.1214/aoms/1177704172
#' @seealso [pxly()]
#' @export
mxmy <- function(x, y, est = median,...){
  out <- est(outer(x, y, FUN = "-"), ...)
  out
}




