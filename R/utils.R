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
  
  # Apply conditional logic for handling values >= 0.5
  data <- data %>%
    mutate(
      normalized_size = ifelse(normalized_size >= 0.5, normalized_size * 1.5, normalized_size)
    )
  
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
generate_alluvium_polygons <- function(flow_data, strata_positions, curve_range = 0.1, n = 50) {
  # Step 1: Check input data structure
  if (!all(c("OmicLayer_from", "stratum_from", "OmicLayer_to", "stratum_to", "size") %in% colnames(flow_data))) {
    stop("Input 'flow_data' is missing required columns.")
  }
  if (!all(c("Layer", "group", "xmin", "xmax", "ymin", "ymax") %in% colnames(strata_positions))) {
    stop("Input 'strata_positions' is missing required columns.")
  }

  # Step 2: Apply layout function
  flow_data_positioned <- layout_flows_within_strata(flow_data, strata_positions)

  # Step 3: Prepare data for polygon generation
  polygon_prep_data <- flow_data_positioned %>%
    select(
      start_x,
      end_x,
      OmicLayer_from,
      stratum_from,
      start_ymin,
      start_ymax,
      OmicLayer_to,
      end_ymin,
      end_ymax,
      stratum_to,
      Drug
    ) %>%
    rename(x = OmicLayer_from) %>%
    mutate(group_id = interaction(stratum_from, x, stratum_to, OmicLayer_to, Drug),fill=Drug) # Create a unique group ID

  # Function to generate polygon points
  generate_polygon <- function(row, n, curve_range) {
    x_start_numeric <- row$start_x 
    x_end_numeric <- row$end_x
    y_start_min <- row$start_ymin
    y_start_max <- row$start_ymax
    y_end_min <- row$end_ymin
    y_end_max <- row$end_ymax

    curve_control_start_y_top <- y_start_max + curve_range
    curve_control_end_y_top <- y_end_max + curve_range
    curve_control_start_y_bottom <- y_start_min - curve_range
    curve_control_end_y_bottom <- y_end_min - curve_range

    # Generate top curve points
    path_top <- bezier::bezier(
      t = seq(0, 1, length.out = n),
      p = matrix(c(x_start_numeric, y_start_max,
                   x_start_numeric + 0.5, curve_control_start_y_top,
                   x_end_numeric - 0.5, curve_control_end_y_top,
                   x_end_numeric, y_end_max),
                 ncol = 2, byrow = TRUE)
    ) %>% as.data.frame() %>% setNames(c("x", "y")) %>% mutate(side = "top")

    # Generate bottom curve points
    path_bottom <- bezier::bezier(
      t = seq(0, 1, length.out = n),
      p = matrix(c(x_start_numeric, y_start_min,
                   x_start_numeric + 0.5, curve_control_start_y_bottom,
                   x_end_numeric - 0.5, curve_control_end_y_bottom,
                   x_end_numeric, y_end_min),
                 ncol = 2, byrow = TRUE)
    ) %>% as.data.frame() %>% setNames(c("x", "y")) %>% mutate(side = "bottom")

    # Combine top and reversed bottom path to form the polygon outline
    polygon_path <- rbind(
      path_top,
      path_bottom %>% arrange(desc(x))
    ) %>%
      mutate(group = row$group_id, fill = row$fill)

    return(polygon_path) %>% select(x, y, group, fill)
  }
  polygon_data <- polygon_prep_data %>%
    rowwise() %>%
    do(generate_polygon(., n = n, curve_range = curve_range)) %>%
    ungroup()

  return(polygon_data)
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
  # Debugging: check the incoming data
  # Step 1: Create two lists of subdataframes based on the conditions
  # List 1: Subdataframes grouped by Drug, OmicLayer_from, and stratum_from
  list_1 <- flow_data %>%
    group_by(Drug, OmicLayer_from, stratum_from) %>%
    group_split()


  # List 2: Subdataframes grouped by Drug, OmicLayer_to, and stratum_to
  list_2 <- flow_data %>%
    group_by(Drug, OmicLayer_to, stratum_to) %>%
    group_split()

  # Initialize an empty list to store the processed sub-dataframes from list_1
  processed_list_1 <- list()
  for (df in list_1) {
    if (nrow(df) > 0) {
      df_processed <- df %>%
        mutate(
          start_x = sapply(OmicLayer_from, function(layer) {
            strata_positions %>%
              filter(Layer == layer) %>%
              pull(xmax) %>%
              first() # Ensure we get a single value
          }),
          end_x = sapply(OmicLayer_to, function(layer) {
            strata_positions %>%
              filter(Layer == layer) %>%
              pull(xmin) %>%
              first() # Ensure we get a single value
          }),

          start_ymin = sapply(1:length(.), function(i) {
            strata_positions[which(.[["OmicLayer_from"]][i] == strata_positions[["Layer"]] & .[["stratum_from"]][i] == strata_positions[["group"]]),] %>%
              pull(ymin) %>%
              first() # Ensure we get a single value
          }) %>% first(),
          start_ymax = sapply(1:length(.), function(i) {
            strata_positions[which(.[["OmicLayer_from"]][i] == strata_positions[["Layer"]] & .[["stratum_from"]][i] == strata_positions[["group"]]),] %>%
              pull(ymax) %>%
              first() # Ensure we get a single value
          }) %>% first(),
               total_size = sum(size)


        ) %>% as.data.frame()
      processed_list_1[[length(processed_list_1) + 1]] <- df_processed
    }
  }
  # Combine the processed sub-dataframes from list_1 back into one
  processed_flow_data_1 <- bind_rows(processed_list_1)

  # Initialize an empty list to store the processed sub-dataframes from list_2
  processed_list_2 <- list()
  for (df in list_2) {
    if (nrow(df) > 0) {
      df_processed <- df %>%
        mutate(
          end_ymin = sapply(1:length(.), function(i) {
            strata_positions[which(.[["OmicLayer_to"]][i] == strata_positions[["Layer"]] & .[["stratum_to"]][i] == strata_positions[["group"]]),] %>%
              pull(ymin) %>%
              first() # Ensure we get a single value
          }) %>% first(),
          end_ymax = sapply(1:length(.), function(i) {
            strata_positions[which(.[["OmicLayer_to"]][i] == strata_positions[["Layer"]] & .[["stratum_to"]][i] == strata_positions[["group"]]),] %>%
              pull(ymax) %>%
              first() # Ensure we get a single value
          }) %>% first(),



        )
      processed_list_2[[length(processed_list_2) + 1]] <- df_processed
    } else {
      processed_list_2[[length(processed_list_2) + 1]] <- df # Keep empty dataframes if any
    }
  }
  # Combine the processed sub-dataframes from list_2 back into one
  processed_flow_data_2 <- bind_rows(processed_list_2)

  # It's important to merge or combine the information from both processed lists
  # based on a common key (e.g., all columns except the newly created ones).
  # A simple way, assuming the order is maintained (which might not be reliable
  # after group_split), is to add the columns from processed_flow_data_2 to
  # processed_flow_data_1 based on matching rows. A safer way would be to perform a join.

  # Let's perform a left_join based on all common columns
  common_cols <- intersect(names(processed_flow_data_1), names(processed_flow_data_2))
  flow_data_updated <- processed_flow_data_1 %>%
    left_join(processed_flow_data_2, by = common_cols)

  # Debugging: check the structure of the updated flow_data
  return(flow_data_updated)
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

  ribbon_data <- generate_alluvium_polygons(flow_data, strata_positions, curve_range = 0.0005)
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
