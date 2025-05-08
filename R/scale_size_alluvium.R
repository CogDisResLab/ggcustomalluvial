#' Scale for Alluvium Size
#'
#' Custom scale function for alluvium size, used in `geom_alluvium_custom()`.
#' It scales the `size` aesthetic according to a specified range.
#'
#' @param range A numeric vector of length 2 specifying the min and max size of the alluvium flows.
#' @param ... Additional arguments passed to `ggplot2`'s scale functions.
#' @return A scale for the alluvium size aesthetic.
#' @export
scale_size_alluvium <- function(range = c(1, 10), ...) {
  ggplot2::scale_size_continuous(range = range, ...)
}
