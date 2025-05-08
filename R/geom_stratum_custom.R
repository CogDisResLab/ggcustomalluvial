#' Custom Geom for Strata
#'
#' A custom geom to represent strata as rectangles. This allows the user to
#' define the positions of strata using `xmin`, `xmax`, `ymin`, and `ymax`.
#'
#' @param mapping Aesthetic mappings (e.g., `aes(xmin, xmax, ymin, ymax)`).
#' @param data A data frame.
#' @param ... Other parameters passed to `ggplot2::geom_rect()`.
#' @return A ggplot layer to visualize strata.
#' @export
geom_stratum_custom <- function(mapping = NULL, data = NULL, ...) {
  ggplot2::geom_rect(mapping = mapping, data = data, ...)
}
