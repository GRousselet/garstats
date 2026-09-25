#' Bus time delays in minutes
#'
#' Bus time delays in minutes for the X10 line from Strathblane to Glasgow.
#' Delays are rounded to the nearest minute. I check my watch when the bus doors open.
#' Positive values mean that the bus was late by x minutes. Negative values indicate that the bus was early.
#' Morning and afternoon bus times are saved separately and do not match, as sometimes I take the bus only
#' one way. To keep things simple i did not record any other details, such as year, time of year, road work, weather.
#'
#' @format ## `bustimes`
#' A data frame with 204 rows and 2 columns:
#' \describe{
#'  \item{period}{Morning or afternoon}
#'  \item{delay}{Rounded delays in minutes. Negative means early, positive means late.}
#'  }
#' @source Data collected by Guillaume A. Rousselet on his occasional bus journeys to work.
#' @name bustimes
#' @docType data
NULL
