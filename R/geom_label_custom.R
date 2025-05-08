#' Custom Geom for Labels
#'
#' A custom geom to represent text labels in a plot. It allows users to add labels
#' with a background box (optional) and customizable appearance.
#'
#' @param mapping Aesthetic mappings for text position and label aesthetics (e.g., `aes(x, y, label)`).
#' @param data A data frame.
#' @param ... Additional parameters passed to `ggplot2::geom_label()`.
#' @return A ggplot layer to visualize custom labels.
#' @export
geom_label_custom <- function(mapping = NULL, data = NULL, ...) {
  ggplot2::geom_label(mapping = mapping, data = data, ...)
}
