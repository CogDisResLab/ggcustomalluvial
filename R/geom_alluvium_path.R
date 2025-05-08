#' Custom Geom for Alluvium Path (Curved Ribbons)
#'
#' A custom geom to represent alluvium flows with curves between strata using polygons.
#' Requires precomputed x, y coordinates forming a smooth path (e.g., from bezier interpolation).
#'
#' @param mapping Aesthetic mappings (e.g., `aes(x, y, group, fill)`).
#' @param data A data frame containing polygon points for the ribbon.
#' @param ... Additional parameters passed to `ggplot2::geom_polygon()`.
#' @return A ggplot layer to visualize curved alluvium ribbons.
#' @export
geom_alluvium_path <- function(mapping = NULL, data = NULL, ...) {
  ggplot2::geom_polygon(mapping = mapping, data = data, ...)
}
