#' Layout Strata Positions
#'
#' This function computes the positions for strata based on a grouping variable.
#' It can automatically distribute strata along the x-axis and calculate y-positions based on the data.
#'
#' @param data A data frame containing the strata data.
#' @param group_var The grouping variable used to define strata.
#' @return A data frame with the computed strata positions.
#' @export
layout_strata_positions <- function(data, group_var) {
  # Calculate the number of strata per layer
  strata_per_layer <- length(unique(data[[group_var]]))  # Number of unique signaling types

  # Create empty columns for the positions
  data$ymin <- NA
  data$ymax <- NA
  data$xmin <- NA
  data$xmax <- NA
  data$Layer_Pos <- NA  # To hold the x positions for layers
  
  # Make sure to retain the 'Signaling' variable
  signaling_column <- data[[group_var]]  # Keep the 'Signaling' column
  
  # Calculate ymin and ymax for each stratum by layer
  # This ensures strata are spaced equally on the y-axis
  for (i in seq_along(unique(data$Layer))) {
    # Get the subset of data corresponding to the current layer
    layer_data <- subset(data, Layer == unique(data$Layer)[i])
    
    # Assign unique y positions for each group within the layer
    n_groups <- nrow(layer_data)
    y_positions <- seq(0, 1, length.out = n_groups)
    
    # Set the ymin and ymax for each stratum
    data$ymin[data$Layer == unique(data$Layer)[i]] <- rev(y_positions)
    data$ymax[data$Layer == unique(data$Layer)[i]] <- rev(y_positions + 0.2) # Padding between strata

    # Calculate xmin and xmax (positions on the x-axis)
    # The xmin and xmax will be determined based on the layer's position on the x-axis
    data$xmin[data$Layer == unique(data$Layer)[i]] <- ((i-1)*10)+4 - 3  # Slight padding to make rectangles visible
    data$xmax[data$Layer == unique(data$Layer)[i]] <- ((i-1)*10)+4 + 3  # Slight padding to make rectangles visible
    
    # Set the layer positions (Layer_Pos) for the strata
    data$Layer_Pos[data$Layer == unique(data$Layer)[i]] <- ((i-1)*10)+4
  }

  # Add the 'Signaling' variable back to the data frame
  data[[group_var]] <- signaling_column

  # Return the data with all necessary variables: xmin, xmax, ymin, ymax, Layer_Pos, and Signaling
  return(data)
}

#' Scale Color for Strata
#'
#' This function scales the color of strata based on a given variable, either continuous or categorical.
#'
#' @param data A data frame containing the strata data.
#' @param group_var The variable to map to color.
#' @param palette The color palette to use (for categorical).
#' @return A vector of colors corresponding to the strata.
#' @export
scale_color_stratum <- function(data, group_var, palette = "Set1") {
  # If the group_var is categorical, use discrete color scale
  if (is.factor(data[[group_var]]) || is.character(data[[group_var]])) {
    library(RColorBrewer)
    color_scale <- scale_color_manual(values = brewer.pal(length(unique(data[[group_var]])), palette))
  } else {
    # For continuous variables, use a gradient scale
    color_scale <- scale_color_gradient(low = "blue", high = "red")
  }
  
  return(color_scale)
}

#' Layout Alluvium Paths
#'
#' This function computes paths for alluvium flows, connecting strata based on some flow data.
#' It calculates the necessary coordinates to draw alluvium paths between strata.
#'
#' @param data A data frame with columns for source and target strata, and the flow value.
#' @param group_var The grouping variable used for the flow paths.
#' @return A data frame with computed x/y positions for alluvium paths.
#' @export
layout_alluvium_paths <- function(data, group_var) {
  # Placeholder: You could define the specific logic for computing the flow paths here
  # For simplicity, we'll create straight paths from one stratum to the next.
  
  # Generate simple x/y paths based on strata positions
  strata_positions <- layout_strata_positions(data, group_var)
  
  # Placeholder: simple straight paths between strata
  data$start_x <- strata_positions$xmin[match(data$source, strata_positions[[group_var]])]
  data$end_x <- strata_positions$xmax[match(data$target, strata_positions[[group_var]])]
  data$start_y <- strata_positions$ymin[match(data$source, strata_positions[[group_var]])]
  data$end_y <- strata_positions$ymax[match(data$target, strata_positions[[group_var]])]
  
  return(data)
}

#' Create Strata Data
#'
#' This function generates a simple data frame for strata, useful for testing or quick visualization.
#'
#' @param n_strata The number of strata to create.
#' @param group_var The variable to use for grouping the strata.
#' @param grouping_levels A vector of values for the grouping variable.
#' @return A data frame with strata positions and grouping variable.
#' @export
create_strata_data <- function(n_axes, group_var = "group", grouping_levels = NULL) {
  if (is.null(grouping_levels)) {
    grouping_levels <- paste("Group", 1:n_strata)
  }
  
  # Create a data frame with strata positions and grouping variable
  data <- data.frame(
    group = rep(grouping_levels, 3),
    xmin = rep(seq(1, n_axes, by = 2), each = n_axes),
    xmax = rep(seq(2, n_axes, by = 2), each = n_axes),
    Layer = c(rep(4, length(grouping_levels)), rep(10, length(grouping_levels)), rep(16, length(grouping_levels)))
  )
  
  return(data)
}

#' Normalize Alluvium Sizes
#'
#' This function normalizes the sizes of alluvium flows to a specific range, typically used for scaling flow width.
#'
#' @param data A data frame with flow data.
#' @param size_var The column in the data frame that represents the size of each flow.
#' @param range The range to which the flow sizes will be normalized (default is 0 to 1).
#' @return A data frame with normalized flow sizes.
#' @export
normalize_alluvium_sizes <- function(data, size_var, range = c(0, 1), group_vars = NULL) {
  if (!is.null(group_vars)) {
    data <- data %>%
      group_by(across(all_of(group_vars))) %>%
      mutate(
        normalized_size = (get(size_var) - min(get(size_var), na.rm = TRUE)) /
                          (max(get(size_var), na.rm = TRUE) - min(get(size_var), na.rm = TRUE)),
        normalized_size = normalized_size * diff(range) + range[1]
      ) %>%
      ungroup()
  } else {
    # Global normalization (current behavior)
    min_size <- min(data[[size_var]], na.rm = TRUE)
    max_size <- max(data[[size_var]], na.rm = TRUE)
    
    data$normalized_size <- (data[[size_var]] - min_size) / (max_size - min_size)
    data$normalized_size <- data$normalized_size * (diff(range)) + range[1]
  }
  
  return(data)
}

generate_alluvium_ribbons <- function(data, id_col = "group", x_col = "x", y_col = "y", size_col = "size") {
  data <- data[order(data[[x_col]], data[[y_col]]), ]
  data$ymin <- data[[y_col]]
  data$ymax <- data$ymin + data[[size_col]]

  # Return in proper order and add group
  data$group <- data[[id_col]]
  return(data)
}

#' Generate Polygon Data for Curved Alluvium Ribbons (with visible curvature)
#'
#' @param data A data frame with x, xend, y, yend, size, group
#' @param n Number of interpolation points per edge
#' @param curve_strength How much to curve the path (positive number)
#' @return A data frame suitable for `geom_polygon()`
#' Generate Polygon Data for Curved Alluvium Ribbons (Cubic Bézier)
#'
#' @param data A data frame with x, xend, y, yend, size, group
#' @param n Number of interpolation points per edge
#' @param curve_strength Vertical offset for control points
#' @return A data frame suitable for `geom_polygon()`
generate_alluvium_polygons <- function(data, n = 50, curve_range = 10) {
  make_sigmoid <- function(range) {
    function(t) {
      s <- 1 / (1 + exp(-range * (t - 0.5)))
      (s - min(s)) / (max(s) - min(s))
    }
  }

  curve_fun <- make_sigmoid(curve_range)
  t_vals <- seq(0, 1, length.out = n)
  f_vals <- curve_fun(t_vals)

  polygons <- lapply(seq_len(nrow(data)), function(i) {
    row <- data[i, ]

    x_path <- row$start_x + (row$end_x - row$start_x) * t_vals
    y_btm <- row$start_y + (row$end_y - row$start_y) * f_vals
    y_top <- y_btm + row$flow_width

    data.frame(
      x = c(x_path, rev(x_path)),
      y = c(y_btm, rev(y_top)),
      group = paste0("ribbon_", i),
      fill = row$Drug
    )
  })

  do.call(rbind, polygons)
}
#' Normalize Flow Sizes
#'
#' This function normalizes the flow sizes to a specific range, useful for adjusting the width of the flow lines.
#'
#' @param data A data frame with flow data.
#' @param size_var The column in the data frame that represents the size of each flow.
#' @param range The range to which the flow sizes will be normalized (default is 0 to 1).
#' @return A data frame with normalized flow sizes.
normalize_flow_sizes <- function(data, size_var, range = c(0, 1)) {
  # Normalize the size values to the desired range
  min_size <- min(data[[size_var]], na.rm = TRUE)
  max_size <- max(data[[size_var]], na.rm = TRUE)
  
  data$normalized_size <- (data[[size_var]] - min_size) / (max_size - min_size)
  data$normalized_size <- data$normalized_size * (diff(range)) + range[1]  # Apply range scaling
  
  return(data)
}

#' Layout Flow Paths
#'
#' This function computes the coordinates for flow paths based on the data, connecting strata
#' based on the flow information.
#'
#' @param data A data frame containing flow information.
#' @param group_var The grouping variable for strata.
#' @return A data frame with computed x/y positions for the flow paths.
layout_flows_within_strata <- function(flow_data, strata_positions) {
  library(dplyr)

  # Join positions for source strata
  flow_data <- flow_data %>%
    left_join(
      strata_positions,
      by = c("OmicLayer_from" = "Layer", "stratum_from" = "group")
    ) %>%
    rename(
      start_x = xmax,
      start_ymin = ymin,
      start_ymax = ymax
    )

  # Prepare target strata with renamed position columns
  strata_positions_target <- strata_positions %>%
    rename(
      end_x = xmin,
      end_ymin = ymin,
      end_ymax = ymax
    )

  # Join positions for target strata
  flow_data <- flow_data %>%
    left_join(
      strata_positions_target,
      by = c("OmicLayer_to" = "Layer", "stratum_to" = "group")
    )

  # Filter out zero-size flows BEFORE allocation
  flow_data <- flow_data %>%
    filter(size > 0)

  # Allocate flows within source strata (from top to bottom),
  # grouped by source layer + stratum only
  flow_data <- flow_data %>%
    group_by(OmicLayer_from, stratum_from) %>%
    arrange(Drug, OmicLayer_to, stratum_to) %>%
    mutate(
      total_start_height = start_ymax - start_ymin,
      n_flows = n(),
      flow_height = total_start_height / n_flows,
      flow_index = row_number(),
      flow_start_y = start_ymax - flow_height * flow_index
    ) %>%
    ungroup()

  # Allocate flows within target strata (from top to bottom),
  # grouped by target layer + stratum only
  flow_data <- flow_data %>%
    group_by(OmicLayer_to, stratum_to) %>%
    arrange(Drug, OmicLayer_from, stratum_from) %>%
    mutate(
      total_end_height = end_ymax - end_ymin,
      n_flows = n(),
      flow_height = total_end_height / n_flows,
      flow_index = row_number(),
      flow_end_y = end_ymax - flow_height * flow_index
    ) %>%
    ungroup()

  # Assign final coordinates
  flow_data$start_y <- flow_data$flow_start_y
  flow_data$end_y <- flow_data$flow_end_y
  flow_data$flow_width <- flow_data$flow_height

  return(flow_data)
}
#' Generate Flow Visualization Data
#'
#' This function prepares the data for visualization by calculating the normalized flow widths.
#' It also determines the path coordinates for each flow.
#'
#' @param data A data frame containing flow information.
#' @param size_col The column representing the flow size.
#' @return A data frame with flow paths and normalized sizes for visualization.
generate_flow_viz_data <- function(data, size_col = "size") {
  # Preserve relevant grouping variable like Drug
  if (!"Drug" %in% colnames(data)) stop("Input data must contain 'Drug' column")

  # Normalize flow sizes
  data <- normalize_flow_sizes(data, size_col)

  # Layout flow paths
  data <- layout_flow_paths(data, "OmicLayer_from")

  # Generate additional columns for visualization
  data$flow_width <- data$normalized_size * 5
  data$flow_y_center <- (data$start_y + data$end_y) / 2

  return(data)
}

#' Plot Flow Paths
#'
#' This function plots the flow paths with the given data, mapping the flow size to the width of the paths.
#'
#' @param data A data frame with flow path and visualization data.
#' @param x_col The column containing x-axis positions (start and end).
#' @param y_col The column containing y-axis positions (start and end).
#' @param width_col The column representing flow width (normalized size).
#' @return A ggplot object with the flow visualization.
plot_flow_paths <- function(data, x_col = "start_x", y_col = "flow_y_center", width_col = "flow_width") {
  library(ggplot2)
  
  ggplot(data, aes(x = !!sym(x_col), y = !!sym(y_col), group = interaction(OmicLayer_from, OmicLayer_to))) +
    geom_segment(aes(xend = end_x, yend = flow_y_center, size = !!sym(width_col)), lineend = "round") +
    scale_size_continuous(range = c(0.1, 5)) +  # Adjust size scale for clarity
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Omic Layer", y = "Flow Path") +
    facet_wrap(~Drug)
}

create_background_strata <- function(n_axes = 3) {
  # Define signaling groups
  signaling_groups <- c(
    "Dopaminergic",
    "Serotonergic",
    "Adrenergic",
    "Glutamatergic",
    "Gabaergic",
    "Non-Canonical Signaling"
  )

  # Step 1: Create strata data across layers
  strata_data <- create_strata_data(
    n_axes = n_axes,
    group_var = "Signaling",
    grouping_levels = signaling_groups
  )

  # Step 2: Layout strata positions
  strata_layout <- layout_strata_positions(strata_data, group_var = "group")

  # Step 3: Plot background strata outlines with labels
  p <- ggplot(strata_layout) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "white", color = "black", linewidth = 0.3
    ) +
    geom_text(
      aes(x = Layer_Pos, y = (ymin + ymax) / 2, label = group),
      size = 3, hjust = 0.5
    ) +
    theme_void() +
    theme(legend.position = "none")

  return(p)
}

generate_alluvium_polygons <- function(data, n = 50, curve_range = 10) {
  make_sigmoid <- function(range) {
    function(t) {
      s <- 1 / (1 + exp(-range * (t - 0.5)))
      (s - min(s)) / (max(s) - min(s))
    }
  }

  curve_fun <- make_sigmoid(curve_range)
  t_vals <- seq(0, 1, length.out = n)
  f_vals <- curve_fun(t_vals)

  polygons <- lapply(seq_len(nrow(data)), function(i) {
    row <- data[i, ]

    x_path <- row$start_x + (row$end_x - row$start_x) * t_vals
    y_btm <- row$start_y + (row$end_y - row$start_y) * f_vals
    y_top <- y_btm + row$flow_width

    data.frame(
      x = c(x_path, rev(x_path)),
      y = c(y_btm, rev(y_top)),
      group = paste0("ribbon_", i),
      fill = row$Drug
    )
  })

  do.call(rbind, polygons)
}

plot_alluvial_from_data <- function(input_data) {
  library(dplyr)
  library(ggplot2)

  # STEP 1: Extract strata positions
  strata_from <- input_data %>%
    select(Layer = OmicLayer_from, group = stratum_from) %>%
    distinct()
  input_data[, c('xstart', 'xstop')] <- NA
  input_data[which(input_data$Layer=="Proteomic"),]$start_x <- input_data[which(input_data$Layer=="Proteomic"),]$xmax
  input_data[which(input_data$Layer=="Proteomic"),]$stop_x <- input_data[which(input_data$Layer=="Proteomic"),]$xmax + 4
  strata_to <- input_data %>%
    select(Layer = OmicLayer_to, group = stratum_to) %>%
    distinct()

  all_strata <- bind_rows(strata_from, strata_to) %>% distinct()
  strata_positions <- layout_strata_positions(all_strata, group_var = "group")

  # STEP 2: Layout flow positions within each stratum
  flow_data <- layout_flows_within_strata(input_data, strata_positions)

  # STEP 3: Construct group and visual properties
  flow_data <- flow_data %>%
    mutate(
      group = paste(Drug, stratum_from, stratum_to, sep = "_"),
      fill = Drug
    )

  # STEP 4: Generate polygon ribbons
  ribbon_data <- generate_alluvium_polygons(flow_data, curve_range = 5)

  # STEP 5: Create background
  background_plot <- create_background_strata(n_axes = 3)

  # STEP 6: Combine plot
  final_plot <- background_plot +
    geom_polygon(
      data = ribbon_data,
      aes(x = x, y = y, group = group, fill = fill),
      color = NA,
      alpha = 0.6
    ) +
    scale_fill_manual(values = c("AH236" = "#1f77b4", "IS141" = "#ff7f0e")) +
    theme_void()

  return(final_plot)
}
