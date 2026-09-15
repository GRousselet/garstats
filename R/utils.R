# Nice colour palette
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cpalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#' Package-specific ggplot2 theme
#'
#' @param axis.title.size Font size for axis titles.
#' @param axis.text.size Font size for axis text.
#' @param plot.title.size Font size for plot titles.
#' @param legend.text.size Font size for legend text.
#' @param legend.title.size Font size for legend titles.
#' @param strip.text.size Font size for facet strip text.
#' @param ... Additional arguments passed to [ggplot2::theme()].
#' @return A complete ggplot2 theme.
#' @export
theme_gar <- function(axis.title.size = 16,
                      axis.text.size = 14,
                      plot.title.size = 20,
                      legend.text.size = 16,
                      legend.title.size = 18,
                      strip.text.size = 20,
                      ...){
  ggplot2::`%+replace%`(ggplot2::theme_bw(), ggplot2::theme(
    plot.title = ggplot2::element_text(size = plot.title.size),
    axis.title = ggplot2::element_text(size = axis.title.size, face = "bold"),
    axis.text = ggplot2::element_text(size = axis.text.size, colour = "black"),
    legend.key.width = grid::unit(1.5, "cm"),
    legend.text = ggplot2::element_text(size = legend.text.size),
    legend.title = ggplot2::element_text(size = legend.title.size),
    strip.text = ggplot2::element_text(size = strip.text.size, face = "bold", colour = "white"),
    strip.background = ggplot2::element_rect(colour = "black", fill = "grey40"),
    complete = TRUE
    ))
}

#' Preserve the order of the names in a variable when converting to a factor
#'
#' @param x A vector to convert to a factor.
#' @return A factor whose levels follow the order of first appearance in `x`.
#' @examples
#' # example code
#' df <- tibble::tibble(cond = c("A", "C", "B"))
#' df$cond <- keeporder(df$cond)
#' @export
keeporder <- function(x){
  x <- as.character(x)
  x <- factor(x, levels=unique(x))
  x
}

#' Print simulation progress
#'
#' @param S Current simulation iteration.
#' @param nsim Total number of simulation iterations.
#' @param inc Frequency at which to print progress updates.
#' @param name Label printed with the initial progress message.
#' @return Write progress to the console. No object is returned.
#' @examples
#' n_iter <- 2000 # simulation iterations
#' inc.step <- 100 # increments
#' for (iter in 1:n_iter) {
#'   sim_counter(iter, n_iter, inc = inc.step, name = "Normal")
#' }
#' @export
sim_counter <- function(S, nsim, inc, name = "Simulation"){
  if(S == 1){
    # print(paste(nsim,"iterations:",S))
    cat(name, nsim,"iterations:",S)
  }
  if(S %% inc == 0){
    # print(paste("iteration",S,"/",nsim))
    cat(" /",S)
  }
}


