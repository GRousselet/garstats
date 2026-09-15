# Quantile-based measures of location, spread, skewness, and kurtosis, using the Harrell-Davis quantile estimator.
# Reference:
# Prendergast, L. A., Dedduwakumara, S., & Staudte, R. G. (2024).
# rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles (arXiv:2410.11093; Version 1).
# arXiv. https://doi.org/10.48550/arXiv.2410.11093

#' Measure of location (default central tendency) relative to a bound (upper or lower), normalised over a reference interval.
#'
#' Implement the distributional shift measure of Locey & Stein (2025), using the Harrell-Davis quantile estimator.
#'
#' @param x Numeric vector.
#' @param q.C Quantile to use as a measure of location. Default 0.5.
#' @param q.L Quantile to use as a lower bound. Default 0.1. Ignored if L and U are supplied.
#' @param q.U Quantile to use as a upper bound. Default 0.9. Ignored if L and U are supplied.
#' @param L Optional fixed lower bound, used instead of q.L. Must be supplied together with U.
#' @param U Optional fixed upper bound, used instead of q.U. Must be supplied together with L.
#' @return Shift measure in the range 0-1.
#' @details
#' The quantity S = (U-C) / (U-L) is a proportion `[0, 1]`.
#' It is related to other quantile measures of skewness in some situations.
#' L = lower bound; U = upper bound.
#' `[L, U]` is a closed interval, with L < U. For observed data the range could be used, but it could be a smaller or larger reference interval.
#' C = measure of location, typically central tendency. Here the default is the Harrell-Davis estimate of the 0.5 quantile (Q.5).
#' By default, S = (Q.9 - Q.5) / (Q.9 - Q.1).
#' S could be re-expressed over `[-1, 1]` by the transformation 2(S-0.5).
#' @references Locey, K. J., & Stein, B. D. (2025).
#' Measurement and Comparison of Distributional Shift with Applications to Ecology, Economics, and Image Analysis.
#' Journal of Statistical Theory and Practice, 19(4), 69. https://doi.org/10.1007/s42519-025-00475-x
#' @export
hd_crub <- function(x, q.C = 0.5, q.L = 0.1, q.U = 0.9, L = NULL, U = NULL){
  if (xor(is.null(L), is.null(U))) stop("L and U must be supplied together")
  if (is.null(L)){
    L <- hd(x, q = q.L)
    U <- hd(x, q = q.U)
  } else if (L >= U) stop("L must be < U")
  C <- hd(x, q = q.C)
  S <- (U-C) / (U-L)
  S
}

#' Lower tail spread
#'
#' Computes the difference between a central and lower quantile.
#'
#' @param x Numeric vector.
#' @param q.C Central quantile. Default = 0.5.
#' @param q.L Lower quantile. Default = 0.1.
#' @return Numeric lower tail spread.
#' @references Prendergast, L. A., Dedduwakumara, S., & Staudte, R. G. (2024).
#' rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles
#' (arXiv:2410.11093; Version 1). arXiv. https://doi.org/10.48550/arXiv.2410.11093
#' @export
hd_lts <- function(x, q.C = 0.5, q.L = 0.1){
  S <- hd(x, q = q.C) - hd(x, q = q.L)
  S
}

#' Upper tail spread
#'
#' Computes the difference between an upper and central quantile.
#'
#' @param x Numeric vector.
#' @param q.C Central quantile. Default = 0.5.
#' @param q.U Upper quantile. Default = 0.9.
#' @return Numeric upper tail spread.
#' @references Prendergast, L. A., Dedduwakumara, S., & Staudte, R. G. (2024).
#' rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles
#' (arXiv:2410.11093; Version 1). arXiv. https://doi.org/10.48550/arXiv.2410.11093
#' @export
hd_uts <- function(x, q.C = 0.5, q.U = 0.9){
  S <- hd(x, q = q.U) - hd(x, q = q.C)
  S
}

#' Interquantile range
#'
#' Computes Q(1 - p) - Q(p).
#'
#' @param x Numeric vector.
#' @param p Lower quantile probability, between 0 and 0.5. Default = 0.25.
#' @return Numeric interquantile range.
#' @references Prendergast, L. A., Dedduwakumara, S., & Staudte, R. G. (2024).
#' rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles
#' (arXiv:2410.11093; Version 1). arXiv. https://doi.org/10.48550/arXiv.2410.11093
#' @export
hd_iqr <- function(x, p = 0.25){
  S <- hd(x, q = 1-p) - hd(x, q = p)
  S
}

#' Normalised interquantile range
#'
#' Computes (Q(1 - p) - Q(p)) / (Q(q.U) - Q(q.L)).
#'
#' @param x Numeric vector.
#' @param p Lower quantile probability, between 0 and 0.5. Default = 0.25.
#' @param q.L Lower reference quantile. Default = 0.1.
#' @param q.U Upper reference quantile. Default = 0.9.
#' @return Numeric normalised interquantile range.
#' @references Prendergast, L. A., Dedduwakumara, S., & Staudte, R. G. (2024).
#' rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles
#' (arXiv:2410.11093; Version 1). arXiv. https://doi.org/10.48550/arXiv.2410.11093
#' @export
hd_niqr <- function(x, p = 0.25, q.L = 0.1, q.U = 0.9){
  S <- (hd(x, q = 1-p) - hd(x, q = p)) / (hd(x, q = q.U) - hd(x, q = q.L))
  S
}

#' Quantile coefficient of dispersion
#'
#' Computes (Q(1 - p) - Q(p)) / (Q(1 - p) + Q(p)).
#'
#' @param x Numeric vector with positive quantiles.
#' @param p Lower quantile probability, between 0 and 0.5. Default = 0.25.
#' @return Numeric quantile coefficient of dispersion.
#' @references Prendergast, L. A., Dedduwakumara, S., & Staudte, R. G. (2024).
#' rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles
#' (arXiv:2410.11093; Version 1). arXiv. https://doi.org/10.48550/arXiv.2410.11093
#' @export
hd_qcd <- function(x, p = 0.25){
  S <- (hd(x, q = 1-p) - hd(x, q = p)) / (hd(x, q = 1-p) + hd(x, q = p))
  S
}

#' Robust coefficient of variation
#'
#' Computes 0.75 * (Q(0.75) - Q(0.25)) / Q(0.5).
#'
#' @param x Numeric vector with a nonzero median.
#' @return Numeric robust coefficient of variation.
#' @references Prendergast, L. A., Dedduwakumara, S., & Staudte, R. G. (2024).
#' rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles
#' (arXiv:2410.11093; Version 1). arXiv. https://doi.org/10.48550/arXiv.2410.11093
#' @export
hd_rcv <- function(x){
  S <- 0.75 * (hd(x, q = 0.75) - hd(x, q = 0.25)) / hd(x, q = 0.5)
  S
}

#' Quantile skewness
#'
#' Computes a generalised quantile skewness measure, with default to Bowley's measure.
#'
#' @param x Numeric vector.
#' @param p Lower quantile probability, between 0 and 0.5. Default = 0.25.
#' @return Numeric quantile skewness.
#' @references Prendergast, L. A., Dedduwakumara, S., & Staudte, R. G. (2024).
#' rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles
#' (arXiv:2410.11093; Version 1). arXiv. https://doi.org/10.48550/arXiv.2410.11093
#' @export
hd_skew <- function(x, p = 0.25){
  S <- (hd(x, q = 1-p) + hd(x, q = p) - 2*hd(x, q = 0.5)) / (hd(x, q = 0.75) - hd(x, q = 0.25))
  S
}

#' Moors quantile kurtosis
#'
#' Computes the quantile-based kurtosis measure of Moors (1988).
#'
#' @param x Numeric vector.
#' @return Numeric Moors quantile kurtosis.
#' @references Prendergast, L. A., Dedduwakumara, S., & Staudte, R. G. (2024).
#' rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles
#' (arXiv:2410.11093; Version 1). arXiv. https://doi.org/10.48550/arXiv.2410.11093
#' @export
hd_kurt <- function(x){
  S <- (hd(x, q = 7/8) - hd(x, q = 5/8) + hd(x, q = 3/8) - hd(x, q = 1/8)) / (hd(x, q = 6/8) - hd(x, q = 2/8))
  S
}

#' Left tail weight
#'
#' Computes a quantile-based measure of left tail weight.
#'
#' @param x Numeric vector.
#' @param p Quantile probability between 0 and 0.5. Default = 0.25.
#' @return Numeric left tail weight.
#' @references Prendergast, L. A., Dedduwakumara, S., & Staudte, R. G. (2024).
#' rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles
#' (arXiv:2410.11093; Version 1). arXiv. https://doi.org/10.48550/arXiv.2410.11093
#' @export
hd_ltw <- function(x, p = 0.25){
  S <- (hd(x, q = (1-p)/2) + hd(x, q = p/2) - 2 * hd(x, q = 0.25)) / (hd(x, q = (1-p)/2) - hd(x, q = p/2))
  S
}

#' Right tail weight
#'
#' Computes a quantile-based measure of right tail weight.
#'
#' @param x Numeric vector.
#' @param p Quantile probability between 0.5 and 1. Default = 0.75.
#' @return Numeric right tail weight.
#' @references Prendergast, L. A., Dedduwakumara, S., & Staudte, R. G. (2024).
#' rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles
#' (arXiv:2410.11093; Version 1). arXiv. https://doi.org/10.48550/arXiv.2410.11093
#' @export
hd_rtw <- function(x, p = 0.75){
  S <- (hd(x, q = (1+p)/2) + hd(x, q = (1-p)/2) - 2 * hd(x, q = 0.75)) / (hd(x, q = (1+p)/2) - hd(x, q = (1-p)/2))
  S
}

#' Upper-to-lower quantile ratio as a simple measure of inequality
#'
#' Computes Q(1 - p) / Q(p), such as the 90/10 ratio.
#'
#' @param x Numeric vector with positive lower quantile.
#' @param p Lower quantile probability, between 0 and 0.5. Default = 0.1.
#' @return Numeric upper-to-lower quantile ratio.
#' @references Prendergast, L. A., Dedduwakumara, S., & Staudte, R. G. (2024).
#' rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles
#' (arXiv:2410.11093; Version 1). arXiv. https://doi.org/10.48550/arXiv.2410.11093
#' @export
hd_ineq <- function(x, p = 0.1){
  S <- hd(x, q = 1-p) / hd(x, q = p)
  S
}
