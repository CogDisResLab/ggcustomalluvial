#' Stat for Custom Stratum Layout
#'
#' A custom statistical transformation to calculate the positions of strata for an alluvial-like plot.
#' This stat computes the `xmin`, `xmax`, `ymin`, and `ymax` positions for each stratum, based on
#' the number of strata and a grouping variable.
#'
#' @param mapping Aesthetic mappings (typically `aes(xmin, xmax, ymin, ymax, fill)`).
#' @param data Data frame.
#' @param n_strata Number of strata to generate.
#' @param group_var The grouping variable to define strata (defaults to `group`).
#' @param ... Other parameters passed to `ggplot2` layers.
#' @return A data frame with the calculated strata positions.
#' @export
stat_stratum_layout <- function(mapping = NULL, data = NULL, n_strata = 5, group_var = "group", ...) {
  # Define the statistical transformation
  ggplot2::layer(
    stat = StatStratumLayout,
    data = data,
    mapping = mapping,
    geom = "rect",  # Strata are represented by rectangles
    position = "identity",
    show.legend = FALSE,
    inherit.aes = TRUE,
    params = list(n_strata = n_strata, group_var = group_var, ...)
  )
}

# Define the StatStratumLayout class
StatStratumLayout <- ggplot2::ggproto(
  "StatStratumLayout", ggplot2::Stat,
  
  # This function is responsible for calculating the strata positions
  compute_layer = function(self, data, params, scales, n_strata = 5, group_var = "group", ...) {
    # Layout strata positions using the utility function defined earlier
    strata_data <- layout_strata_positions(data, group_var)
    
    # Create the strata with necessary aesthetics
    strata_data$xmin <- rep(1:n_strata, each = n_strata)
    strata_data$xmax <- strata_data$xmin + 1  # Define the width of each stratum
    strata_data$ymin <- strata_data$ymin * 10  # Scale to fit the plot's y range (adjust as needed)
    strata_data$ymax <- strata_data$ymin + 0.1  # Define stratum height
    
    # Return the computed strata positions
    return(strata_data)
  }
)
