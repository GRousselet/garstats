utils::globalVariables(c("gp", "x", "y"))

#' Plot a bootstrap distribution
#'
#' Draws a density curve for a bootstrap distribution exported by the
#' functions in `wilcoxmod.R`. The `x` argument may be either the numeric
#' `boot.estimates` vector or the complete result list returned by a bootstrap
#' function.
#' @param x Numeric bootstrap estimates, or a result list containing
#'   `boot.estimates`.
#' @param ci Optional length-two confidence interval. If `x` is a result
#'   list, `ci` is taken from `x$ci` when available. May include `-Inf`
#'   or `Inf` for a one-sided interval; its finite bound is drawn as a
#'   full-height dashed vertical line.
#' @param ci_level Confidence level used to compute `ci` from `x` when `ci`
#'   is not supplied.
#' @param alternative Type of interval when `ci` is computed from `ci_level`:
#'   `"two.sided"` (default), `"less"` (upper bound only), or `"greater"`
#'   (lower bound only).
#' @param estimate Optional point estimate. If `x` is a result list,
#'   `estimate` is taken from `x$estimate` when available.
#' @param hyp Optional null-hypothesis value. If `x` is a result list,
#'   `hyp` is taken from `x$hyp` when available.
#' @param show_hyp Logical; draw the null-hypothesis reference line when a
#'   finite `hyp` value is available.
#' @param label_x_offset Horizontal label offset in data units. By default,
#'   a data-dependent offset of 3% of the bootstrap range is used.
#' @param label_y_offset Vertical label offset in density units. By default,
#'   a data-dependent offset of 3% of the density height is used.
#' @param label_height Height of the label position as a fraction of the
#'   density height. Defaults to 0.1.
#' @param show_ci Logical; draw the confidence interval when `ci` is available.
#' @param show_estimate Logical; draw the point estimate when `estimate` is
#'   available.
#' @param xlab,ylab Axis labels. Default to "Bootstrap estimates" and "Density".
#' @param line_colour,ci_colour,estimate_colour Colours for the density,
#'   confidence interval, and estimate.
#' @param ci_linewidth Line width for the confidence interval.
#' @param label_fill Fill colour used for endpoint labels.
#' @param label_colour Text colour used for endpoint labels.
#' @param label_size Label text size.
#' @param one_sided_linetype Line type for a one-sided confidence bound.
#' @param hyp_linetype Line type for the null-hypothesis reference line.
#' @param theme A ggplot2 theme. Defaults to `garstats::theme_gar()`.
#' @param ... Additional arguments passed to `stats::density()`.
#'
#' @return A ggplot object.
#' @export
plot.boot <- function(x, ci = NULL, ci_level = 0.97,
                      alternative = c("two.sided", "less", "greater"),
                      estimate = NULL, hyp = NULL, show_hyp = TRUE,
                      label_x_offset = NULL, label_y_offset = NULL,
                      label_height = 0.1, show_ci = TRUE,
                      show_estimate = TRUE, xlab = NULL, ylab = NULL,
                      line_colour = "black", ci_colour = "orange",
                      ci_linewidth = 2, estimate_colour = "black",
                      label_fill = "orange", label_colour = "white",
                      label_size = 4, one_sided_linetype = "solid",
                      hyp_linetype = "dashed",
                      theme = NULL, ...) {
  alternative <- match.arg(alternative)
  if (is.list(x)) {
    result <- x
    x <- result$boot.estimates
    if (is.null(estimate)) estimate <- result$estimate
    if (is.null(hyp)) hyp <- result$hyp
  }
  if (!is.numeric(x) || length(x) < 2L || any(!is.finite(x))) {
    stop("x must contain at least two finite numeric estimates.")
  }
  if (is.null(ci)) {
    if (!is.numeric(ci_level) || length(ci_level) != 1L || ci_level <= 0 || ci_level >= 1) {
      stop("ci_level must be a single value strictly between 0 and 1.")
    }
    alpha <- 1 - ci_level
    ci <- switch(alternative,
      two.sided = quantile(x, c(alpha / 2, 1 - alpha / 2), names = FALSE),
      less = c(-Inf, quantile(x, ci_level, names = FALSE)),
      greater = c(quantile(x, alpha, names = FALSE), Inf)
    )
  } else {
    if (!is.numeric(ci) || length(ci) != 2L || anyNA(ci)) {
      stop("ci must be a numeric vector of length two without NA values.")
    }
  }
  if (!is.null(estimate) &&
      (!is.numeric(estimate) || length(estimate) != 1L || !is.finite(estimate))) {
    stop("estimate must be one finite numeric value.")
  }
  if (!is.null(hyp) &&
      (!is.numeric(hyp) || length(hyp) != 1L || !is.finite(hyp))) {
    stop("hyp must be one finite numeric value.")
  }
  density_data <- stats::density(x, ...)
  density_df <- data.frame(x = density_data$x, y = density_data$y)
  x_span <- diff(range(density_df$x))
  y_span <- diff(range(density_df$y))
  if (x_span == 0) x_span <- 1
  if (y_span == 0) y_span <- 1
  if (is.null(label_x_offset)) label_x_offset <- 0.03 * x_span
  if (is.null(label_y_offset)) label_y_offset <- 0.03 * y_span
  if (is.null(xlab)) xlab <- "Bootstrap estimates"
  if (is.null(ylab)) ylab <- "Density"
  if (is.null(theme)) theme <- garstats::theme_gar()

  plot <- ggplot2::ggplot(density_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(linewidth = 1.2, colour = line_colour) +
    ggplot2::labs(x = xlab, y = ylab)
  if (show_ci) {
    lower_inf <- !is.finite(ci[1])
    upper_inf <- !is.finite(ci[2])
    if (!(lower_inf && upper_inf)) {
      label_y <- max(density_df$y) * label_height + label_y_offset
      if (lower_inf || upper_inf) { # one-sided interval
        bound <- if (lower_inf) ci[2] else ci[1]
        bound_label <- if (lower_inf) "Upper bound" else "Lower bound"
        plot <- plot + ggplot2::geom_vline(
          xintercept = bound, linewidth = ci_linewidth,
          linetype = one_sided_linetype, colour = ci_colour) +
          ggplot2::annotate(
          "label", x = bound, y = label_y,
          label = paste(bound_label, round(bound, 2)),
          size = label_size, colour = label_colour, fill = ci_colour,
          fontface = "bold"
        )
      } else { # two-sided interval
        plot <- plot + ggplot2::geom_segment(
          x = ci[1], xend = ci[2], y = 0, yend = 0,
          linewidth = ci_linewidth, lineend = "round",
          colour = ci_colour
        )
        plot <- plot + ggplot2::annotate(
          "label", x = ci[1] + label_x_offset, y = label_y,
          label = round(ci[1], 2),
          size = label_size, colour = label_colour, fill = ci_colour,
          fontface = "bold") +
          ggplot2::annotate("label", x = ci[2] - label_x_offset, y = label_y,
                            label = round(ci[2], 2),
                            size = label_size, colour = label_colour, fill = ci_colour,
                            fontface = "bold")
      }
    }
  }
  if (show_estimate && !is.null(estimate)) {
    plot <- plot + ggplot2::geom_vline(xintercept = estimate,
                                       colour = estimate_colour, linewidth = 0.8)
  }
  if (show_hyp && !is.null(hyp)) {
    plot <- plot + ggplot2::geom_vline(xintercept = hyp,
                                       linetype = hyp_linetype, linewidth = 0.8)
  }
  plot + theme
}

 #' Plot paired bootstrap estimates with a confidence ellipse
 #'
 #' @param x Result list from `twosampb()` or a numeric vector of bootstrap
 #' estimates from group 1.
 #' @param y Numeric bootstrap estimates from group 2, when `x` is numeric.
 #' @param level Confidence level for `ggplot2::stat_ellipse()`.
 #' @param ellipse_linetype Ellipse line type.
 #' @param ellipse_colour Ellipse colour.
 #' @param ellipse_linewidth Ellipse linewidth.
 #' @param plot_diag Default to TRUE to plot identity reference line
 #' @param diag_linetype Identity reference line type.
 #' @param diag_colour Identity reference line colour.
 #' @param diag_linewidth Identity reference line width.
 #' @param point_colour Point colour.
 #' @param point_alpha Point transparency.
 #' @param xlab,ylab Axis labels.
 #' @param plot.theme A ggplot2 theme. Defaults to `garstats::theme_gar()`.
 #' @param ... Additional arguments accepted for compatibility with [plot()].
 #' @return A ggplot object.
 #' @export
plot.boot2 <- function(x, y = NULL, level = 0.97,
                       ellipse_linetype = "dashed", ellipse_colour = "darkorange",
                       ellipse_linewidth = 1,
                       plot_diag = TRUE, diag_linetype = "dashed",
                       diag_colour = "black", diag_linewidth = 0.75,
                       point_colour = "black", point_alpha = 0.35,
                       xlab = "Bootstrap group 1 estimates",
                       ylab = "Bootstrap group 2 estimates",
                       plot.theme = NULL, ...) {
  if (is.list(x)) {
    y <- x$boot.estimates.y
    x <- x$boot.estimates.x
  }
  if (!is.numeric(x) || !is.numeric(y) ||
      length(x) < 3L || length(x) != length(y) ||
      any(!is.finite(x)) || any(!is.finite(y))) {
    stop("x and y must be finite numeric vectors of equal length.")
  }
  if (level <= 0 || level >= 1) stop("level must be between 0 and 1.")
  data <- data.frame(x = x, y = y)
  if (is.null(plot.theme)) plot.theme <- garstats::theme_gar()
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(colour = point_colour, alpha = point_alpha) +
    ggplot2::stat_ellipse(level = level, linetype = ellipse_linetype,
                          colour = ellipse_colour, linewidth = ellipse_linewidth) +
    ggplot2::labs(x = xlab, y = ylab) +
    plot.theme +
    ggplot2::theme(aspect.ratio = 1)
  if (plot_diag){
    plot <- plot +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = diag_colour,
                           linewidth = diag_linewidth, linetype = diag_linetype)
  }
  plot
}

#' Plot one or two empirical cumulative distribution functions
#'
#' @param x Numeric vector for the first function.
#' @param y Optional numeric vector for the second function.
#' @param xlab,ylab Axis labels.
#' @param theme A ggplot2 theme, with default to `garstats::theme_gar()`.
#' @param x.name,y.name Labels for the plotted functions.
#' @param x.col,y.col Colours for the plotted functions.
#' @param x.lt,y.lt Line types for the plotted functions.
#' @param legend.width Width of the legend keys.
#' @param ... Additional arguments accepted for compatibility with [plot()].
#' @return A ggplot2 object.
#' @export
plot.ecdf <- function(x, y = NULL, xlab = "Measurements", ylab = "F(x) = proportion <= x",
                      theme = NULL, x.name = "Group 1", y.name = "Group 2",
                      x.col = "black", y.col = "grey60",
                      x.lt = "longdash", y.lt = "solid",
                      legend.width = 4, ...) {
  if (!is.null(y)) {
    nx <- length(x)
    ny <- length(y)
    df <- tibble::tibble(x = c(x, y),
                 gp = factor(c(rep(x.name, nx), rep(y.name, ny)),
                             levels = c(x.name, y.name)))
    line_types <- c(x.lt, y.lt)
    colours <- c(x.col, y.col)
  } else {
    df <- tibble::tibble(x = x, gp = factor(x.name, levels = x.name))
    line_types <- x.lt
    colours <- x.col
  }
  # make plot
  if (is.null(theme)) theme <- garstats::theme_gar()
  ggplot2::ggplot(df, ggplot2::aes(x=x, colour=gp, linetype=gp)) +
    ggplot2::stat_ecdf(geom = "step", linewidth = 1) +
    ggplot2::labs(x = xlab, y = ylab) +
    theme +
        ggplot2::scale_linetype_manual(values = line_types) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2),
          linetype = ggplot2::guide_legend(nrow = 2)) +
        ggplot2::theme(legend.title = ggplot2::element_blank(),
         legend.key.width = grid::unit(legend.width, "line"),
          legend.position = c(0.02, 0.98),
          legend.justification = c("left", "top"),
          legend.direction = "vertical"
    )
}


