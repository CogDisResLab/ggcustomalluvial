#' Layout Strata Positions with User-Defined Order
#'
#' This function computes the positions for strata based on a user-provided order.
#'
#' @param data A data frame containing the strata data with 'Layer' and 'group' columns.
#' @param strata_order A named list where names are omics layers and values are character vectors
#'                     specifying the desired order of strata within each layer.
#' @return A data frame with the computed strata positions (xmin, xmax, ymin, ymax, Layer_Pos).
#' @export

layout_strata_positions <- function(data, strata_order) {
  library(dplyr)

  # Extract unique layers
  unique_layers <- unique(data$Layer)

  for (layer in unique_layers) {
    # Filter and prepare data for current layer
    layer_data <- data %>%
      filter(Layer == layer) %>%
      mutate(group = trimws(group)) %>%
      select(-any_of(c("ymin", "ymax")))  # Remove potential conflicting columns

    unique_groups <- unique(layer_data$group)

    if (length(unique_groups) > 0) {
      # Determine the order of groups
      ordered_groups <- strata_order[[layer]] %||% rev(sort(unique_groups))

      # Compute vertical layout values
      padding_proportion <- 0.05
      total_padding <- (length(ordered_groups) + 1) * padding_proportion
      available_height <- 1 - total_padding

      if (available_height > 0) {
        stratum_height <- available_height / length(ordered_groups)
        y_positions_start <- seq(padding_proportion, 1 - padding_proportion, length.out = length(ordered_groups)) - (stratum_height / 2)

        group_y_positions <- data.frame(
          group = ordered_groups,
          ymin = y_positions_start,
          ymax = y_positions_start + stratum_height
        )

        # Join new ymin/ymax values into layer_data
        layer_data <- layer_data %>%
          left_join(group_y_positions, by = "group")

        # Merge the updated layer_data back into the full data
        data <- data %>%
          filter(Layer != layer) %>%
          bind_rows(layer_data)
      }
    }

    # Ensure x columns exist before mutate
    if (!"xmin" %in% names(data)) data$xmin <- NA_real_
    if (!"xmax" %in% names(data)) data$xmax <- NA_real_
    if (!"Layer_Pos" %in% names(data)) data$Layer_Pos <- NA_real_

    # Compute x positions
    layer_index <- which(unique_layers == layer)
    data <- data %>%
      mutate(
        xmin = ifelse(Layer == layer, (layer_index - 1) * 10 + 1, xmin),
        xmax = ifelse(Layer == layer, (layer_index - 1) * 10 + 7, xmax),
        Layer_Pos = ifelse(Layer == layer, (layer_index - 1) * 10 + 4, Layer_Pos)
      )
  }
  return(data)
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
  # Extract strata positions
  strata_positions <- layout_strata_positions(data, group_var)

  # Initialize new variables to store the adjusted y positions for start and end
  data$start_y <- strata_positions$ymin[match(data$source, strata_positions[[group_var]])]
  data$end_y <- strata_positions$ymax[match(data$target, strata_positions[[group_var]])]

  # Separate strata positions by source and target group
  source_group_ymin <- strata_positions %>%
    group_by(group) %>%
    summarise(min_y = min(ymin)) %>%
    rename(source_group = group)

  target_group_ymax <- strata_positions %>%
    group_by(group) %>%
    summarise(max_y = max(ymax)) %>%
    rename(target_group = group)

  # Apply the necessary adjustments to avoid overlap
  data <- data %>%
    left_join(source_group_ymin, by = c("source" = "source_group")) %>%
    left_join(target_group_ymax, by = c("target" = "target_group")) %>%
    mutate(
      # Adjust start_y to avoid overlap
      start_y = start_y + 0.02 * (seq_along(start_y) - 1),
      # Adjust end_y similarly
      end_y = end_y - 0.02 * (seq_along(end_y) - 1)
    )

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
create_strata_data <- function(n_axes, grouping_levels) {
  n_groupings <- length(grouping_levels)
  total_rows <- n_axes * n_groupings  # This should be 6 if n_axes = 3 and n_groupings = 2

  # Generate data frame with correct row lengths
  strata_data <- data.frame(
    group = rep(grouping_levels, each = n_axes),
    xmin = rep(seq(1, n_axes, by = 2), length.out = total_rows),
    xmax = rep(seq(2, n_axes + 1, by = 2), length.out = total_rows),
    Layer = rep(4:6, length.out = total_rows))

  return(strata_data)
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

  # Check if required columns exist
  required_cols <- c("start_x", "end_x", "start_ymin", "end_ymin", "start_ymax", "end_ymax", "Drug")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns: ", paste(missing_cols, collapse = ", ")))
  }

  curve_fun <- make_sigmoid(curve_range)
  t_vals <- seq(0, 1, length.out = n)
  f_vals <- curve_fun(t_vals)

  polygons <- lapply(seq_len(nrow(data)), function(i) {
    row <- data[i, ]

    # Check for NA values and handle them
    if (any(is.na(row[c("start_ymin", "end_ymin", "start_ymax", "end_ymax")]))) {
      message(paste("Skipping row", i, "due to missing y positions"))
      return(NULL)  # Skip this iteration if there are missing values
    }

    row$flow_width <- abs(row$start_ymax - row$start_ymin)

    # Use linear X, curved Y (adjusting for the flow width)
    x_path <- row$start_x + (row$end_x - row$start_x) * t_vals
    y_btm  <- row$start_ymin + (row$end_ymin - row$start_ymin) * f_vals
    y_top  <- row$start_ymax + (row$end_ymax - row$start_ymax) * f_vals

    # Create a data frame with the required columns
    data.frame(
      x = c(x_path, rev(x_path)),
      y = c(y_btm, rev(y_top)),
      group = paste0("ribbon_", i),
      fill = factor(rep(row$Drug, length(x_path) * 2), levels = unique(data$Drug))
    )
  })

  # Remove NULL entries caused by skipped rows
  polygons <- polygons[!sapply(polygons, is.null)]

  # Return combined polygons or empty data frame if no valid polygons
  if (length(polygons) > 0) {
    return(do.call(rbind, polygons))
  } else {
    return(data.frame())
  }
}
#' Layout Flow Paths
#'
#' This function computes the coordinates for flow paths based on the data, connecting strata
#' based on the flow information.
#'
#' @param flow_data A data frame containing flow information with columns
#'   like OmicLayer_from, stratum_from, OmicLayer_to, stratum_to, and size.
#' @param strata_positions A data frame with strata positions including
#'   Layer, group, xmin, xmax, ymin, and ymax.
#' @return A data frame with computed x/y positions for the flow paths.
layout_flows_within_strata <- function(flow_data, strata_positions) {
  library(dplyr)


  # Join with starting strata positions
  flow_data <- flow_data %>%
    left_join(strata_positions %>% mutate(Layer = as.character(Layer)),
              by = c("OmicLayer_from" = "Layer", "stratum_from" = "group")) %>%
    rename(start_x = xmax, start_ymin = ymin, start_ymax = ymax) # Alluvia originate from the end of the source stratum


  # Join with ending strata positions - MODIFIED TO CONNECT TO THE BEGINNING (xmin)
  strata_positions_target <- strata_positions %>%
    rename(
      end_x = xmin, # Changed from xmax to xmin
      end_ymin = ymin,
      end_ymax = ymax
    )

  
  flow_data <- flow_data %>%
    left_join(
      strata_positions_target %>% mutate(Layer = as.character(Layer)),
      by = c("OmicLayer_to" = "Layer", "stratum_to" = "group")
    )
  # Ensure zero-size flows are removed
  flow_data <- flow_data %>%
    filter(size > 0)

  return(flow_data)
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

create_background_strata <- function(data, group_var = "group", layer_var = "Layer") {
  strata_data <- data %>%
    select(!!sym(group_var), !!sym(layer_var)) %>%
    distinct() %>%
    rename(group = !!sym(group_var), Layer = !!sym(layer_var))

  strata_layout <- layout_strata_positions(strata_data, group_var = "group")

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

#' Plot Alluvial Diagram from Data with User-Defined Strata Order
#'
#' This function takes flow data and orders for omics layers and strata to
#' generate an alluvial plot.
#'
#' @param input_data A data frame containing flow information.
#' @param omics_order A character vector specifying the order of omics layers.
#' @param strata_order A named list where names are omics layers and values are character vectors
#'                     specifying the desired order of strata within each layer.
#' @return A ggplot object representing the alluvial plot.
#' @export
plot_alluvial_from_data <- function(input_data, omics_order, strata_order = NULL) {
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(tidyr) # Ensure tidyr is loaded

  if (is.null(omics_order)) {
    stop("Error: The 'omics_order' argument must be provided.")
  }

  # Convert omics layers to ordered factors
  input_data <- input_data %>%
    mutate(
      OmicLayer_from = factor(OmicLayer_from, levels = omics_order, ordered = TRUE),
      OmicLayer_to = factor(OmicLayer_to, levels = omics_order, ordered = TRUE)
    )

  # Check for incorrect flow direction
  invalid_flows <- input_data %>%
    filter(as.numeric(OmicLayer_to) <= as.numeric(OmicLayer_from))

  if (nrow(invalid_flows) > 0) {
    invalid_flow_details <- invalid_flows %>%
      select(OmicLayer_from, stratum_from, OmicLayer_to, stratum_to) %>%
      head() %>%
      unite(col = "flow", sep = " -> ") %>%
      pull(flow)

    error_message <- paste(
      "Error: Some flows are moving against or staying within the specified omics layer order.",
      "\nExamples of invalid flows:",
      paste(invalid_flow_details, collapse = ", "),
      if (nrow(invalid_flows) > 5) paste0("... and ", nrow(invalid_flows) - 5, " more.") else ""
    )
    stop(error_message)
  }

  # Collect all strata
  strata_from <- input_data %>%
    select(Layer = OmicLayer_from, group = stratum_from) %>%
    distinct()
  strata_to <- input_data %>%
    select(Layer = OmicLayer_to, group = stratum_to) %>%
    distinct()
  all_strata <- bind_rows(strata_from, strata_to) %>% distinct()

  # Layout strata positions using the provided order
  strata_order_cleaned <- lapply(strata_order, function(x) gsub("\\s+", "", x))

  strata_positions <- layout_strata_positions(all_strata, strata_order)

  flow_data <- layout_flows_within_strata(input_data, strata_positions)

  flow_data <- flow_data %>%
    mutate(
      group = paste(Drug, stratum_from, stratum_to, sep = "_"),
      start_layer_num = as.numeric(factor(OmicLayer_from, levels = omics_order)),
      end_layer_num = as.numeric(factor(OmicLayer_to, levels = omics_order))
    )

  # Defensive renaming: Ensure start_x and end_x don't exist
  if (!("start_x" %in% names(flow_data))) {
    flow_data <- flow_data %>% rename(start_x = start_layer_num)
  } else {
    flow_data <- flow_data %>% rename(start_x_calc = start_layer_num) # Rename to a temp name
  }

  if (!("end_x" %in% names(flow_data))) {
    flow_data <- flow_data %>% rename(end_x = end_layer_num)
  } else {
    flow_data <- flow_data %>% rename(end_x_calc = end_layer_num) # Rename to a temp name
  }

  ribbon_data <- generate_alluvium_polygons(flow_data, curve_range = 0.5)
  unique_drugs <- unique(input_data$Drug)
  ribbon_data$fill <- factor(ribbon_data$fill, levels = unique_drugs)

  # Dynamically assign colors
  drug_colors <- setNames(brewer.pal(length(unique_drugs), "Set2"), unique_drugs)

  final_plot <- ggplot(ribbon_data, aes(x = x, y = y, group = group, fill = fill)) +
    geom_polygon(color = NA, alpha = 0.6) +
    geom_rect(
      data = strata_positions,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "white", color = "black", linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = strata_positions,
      aes(x = Layer_Pos, y = (ymin + ymax) / 2, label = group),
      size = 3, hjust = 0.5,
      inherit.aes = FALSE
    ) +
    scale_x_continuous(breaks = seq_along(omics_order), labels = omics_order) +
    scale_fill_manual(name = "Drug", values = drug_colors, limits = unique_drugs) +
    theme_void() +
    theme(legend.position = "none")

  return(final_plot)
}
