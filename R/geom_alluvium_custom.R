#' Custom Geom for Alluvium (Ribbons)
#'
#' A custom geom to represent alluvium flows as ribbons (like in traditional alluvial plots).
#' Expects aesthetics like `aes(x, ymin, ymax, group, fill)`.
#'
#' @param mapping Aesthetic mappings (e.g., `aes(x, ymin, ymax, group, fill)`).
#' @param data A data frame containing ribbon data.
#' @param ... Additional parameters passed to `ggplot2::geom_ribbon()`.
#' @return A ggplot layer to visualize alluvium flows using ribbons.
#' @export
geom_alluvium_custom <- function(mapping = NULL, data = NULL, ...) {
  ggplot2::geom_ribbon(mapping = mapping, data = data, ...)
}
