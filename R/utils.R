library(purrr)
library(dplyr)
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
      padding_proportion <- 0.00
      total_padding <- (length(ordered_groups) + 1) * padding_proportion
      available_height <- 1 - total_padding

      if (available_height > 0) {
        stratum_height <- available_height / length(ordered_groups)
        y_positions_start <- seq(0, 1 - stratum_height, length.out = length(ordered_groups))

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

#' Normalize a Vector to Proportions
#'
#' This function normalizes a numeric vector so that its values sum to 1.
#' It is designed to be used inside dplyr::mutate() or similar workflows.
#'
#' @param x A numeric vector to normalize.
#'
#' @return A numeric vector where all non-NA values sum to 1.
#' @export
normalize_alluvium_sizes <- function(x) {
  total <- sum(x, na.rm = TRUE)
  if (total == 0 || is.na(total)) {
    warning("Sum of values is zero or NA. Returning NA for all elements.")
    return(rep(NA_real_, length(x)))
  }
  x / total
}

#' Generate Polygon Data for Curved Alluvium Ribbons (Cubic Bézier Curves)
#'
#' @param flow_data Data frame with columns: OmicLayer_from, stratum_from, OmicLayer_to, stratum_to, size, Drug
#' @param strata_positions Data frame with strata layout info: Layer, group, xmin, xmax, ymin, ymax
#' @param curve_range Numeric scalar controlling vertical curve offset (default 0.1)
#' @param n Number of interpolation points per edge (default 50)
#' @return Data frame suitable for geom_polygon(), with columns x, y, group, fill
generate_alluvium_polygons <- function(flow_data, strata_positions, curve_range = 0.1, n = 50) {
  library(dplyr)
  library(bezier)
  flow_pos <- flow_data

  # Validate required columns in input
  required_flow_cols <- c("OmicLayer_from", "stratum_from", "OmicLayer_to", "stratum_to", "size", "Drug")
  if (!all(required_flow_cols %in% names(flow_data))) {
    stop("Input 'flow_data' is missing required columns: ", 
         paste(setdiff(required_flow_cols, names(flow_data)), collapse = ", "))
  }
  required_strata_cols <- c("Layer", "group", "xmin", "xmax", "ymin", "ymax")
  if (!all(required_strata_cols %in% names(strata_positions))) {
    stop("Input 'strata_positions' is missing required columns: ", 
         paste(setdiff(required_strata_cols, names(strata_positions)), collapse = ", "))
  }

  # Prepare data for polygon generation
  polygon_prep <- flow_pos %>%
    mutate(
      group_id = interaction(stratum_from, OmicLayer_from, stratum_to, OmicLayer_to, Drug, drop = TRUE),
      fill = Drug
    ) %>%
    select(
      start_x, end_x,
      start_ymin, start_ymax,
      end_ymin, end_ymax,
      group_id, fill
    )

  # Function to generate polygon points for a single ribbon (row)
  generate_polygon <- function(row, n_points, curve_height) {
    # Top curve points (Bezier curve)
    top_curve <- bezier::bezier(
      t = seq(0, 1, length.out = n_points),
      p = matrix(c(
        row$start_x, row$start_ymax,
        row$start_x + 0.5, row$start_ymax + curve_height,
        row$end_x - 0.5, row$end_ymax + curve_height,
        row$end_x, row$end_ymax
      ), ncol = 2, byrow = TRUE)
    )
    top_df <- as.data.frame(top_curve)
    names(top_df) <- c("x", "y")
    top_df$side <- "top"

    # Bottom curve points (Bezier curve)
    bottom_curve <- bezier::bezier(
      t = seq(0, 1, length.out = n_points),
      p = matrix(c(
        row$start_x, row$start_ymin,
        row$start_x + 0.5, row$start_ymin - curve_height,
        row$end_x - 0.5, row$end_ymin - curve_height,
        row$end_x, row$end_ymin
      ), ncol = 2, byrow = TRUE)
    )
    bottom_df <- as.data.frame(bottom_curve)
    names(bottom_df) <- c("x", "y")
    bottom_df$side <- "bottom"

    # Combine points to form polygon outline (top forward, bottom reversed)
    polygon_df <- bind_rows(
      top_df,
      bottom_df %>% arrange(desc(x))
    )

    polygon_df$group <- row$group_id
    polygon_df$fill <- row$fill

    return(polygon_df %>% select(x, y, group, fill))
  }

  # Generate polygons for all ribbons
  polygon_data <- polygon_prep %>%
    rowwise() %>%
    do(generate_polygon(., n_points = n, curve_height = curve_range)) %>%
    ungroup()

  return(polygon_data)
}

split_positions_by_size <- function(df, output_stratum = FALSE) {
  n <- nrow(df)
  if (n == 0) {
    return(matrix(numeric(0), ncol = 2, nrow = 0))
  }

  vector_list <- matrix(NA_real_, ncol = 2, nrow = n)

  size_col <- "size"
  ymin_col <- ifelse(!output_stratum, "start_ymin", "end_ymin")
  ymax_col <- ifelse(!output_stratum, "start_ymax", "end_ymax")

  grouping_vars <- if (!output_stratum) {
    c("stratum_to", "Drug")
  } else {
    c("stratum_from", "Drug")
  }

  df <- df %>%
    mutate(row_id = row_number()) %>%
    group_by(across(all_of(grouping_vars))) %>%
    arrange(.by_group = TRUE) %>%
    ungroup()

  df_split <- group_split(df, .keep = TRUE)

  for (sub_df in df_split) {
    ids <- sub_df$row_id
    y_max <- sub_df[[ymax_col]][1]
    y_min <- sub_df[[ymin_col]][1]

    stratum_height <- y_max - y_min
    total_size <- sum(sub_df[[size_col]])

    # Defensive: If total_size is zero or tiny, replace all sizes by equal shares
    if (total_size < .Machine$double.eps) {
      rel_sizes <- rep(1 / nrow(sub_df), nrow(sub_df))  # equal shares
    } else {
      rel_sizes <- sub_df[[size_col]] / total_size
    }

    # Now scale relative sizes *exactly* to stratum height
    scaled_heights <- rel_sizes * stratum_height

    current_y <- y_max
    for (i in seq_along(ids)) {
      h <- scaled_heights[i]
      vector_list[ids[i], ] <- c(current_y, current_y - h)
      current_y <- current_y - h
    }
  }

  return(vector_list)
}


apply_split_positions_by_size <- function(df, OmicLayer, stratum, end = FALSE) {
  filtered_df <- if (end) {
    # For end (output_stratum = TRUE), group by stratum_from and Drug, so filter by stratum_from here
    df %>% filter(OmicLayer_to == OmicLayer, stratum_from == stratum)
  } else {
    # For start (output_stratum = FALSE), group by stratum_to and Drug, so filter by stratum_to here
    df %>% filter(OmicLayer_from == OmicLayer, stratum_to == stratum)
  }
  if (nrow(filtered_df) == 0) {
    return(matrix(numeric(0), ncol = 2, nrow = 0))
  }
  split_positions_by_size(filtered_df, output_stratum = end)
}

layout_flows_within_strata <- function(flow_data, strata_positions) {
  # Relative sizes
  flow_data <- flow_data |>
    mutate(relative_size = size / max(size))

  # Join start (from) positions — keep ymin/ymax
  flow_data <- flow_data |>
    left_join(
      strata_positions |>
        select(group, Layer, start_x = xmax, start_ymin = ymin, start_ymax = ymax) |>
        rename(stratum_from = group, OmicLayer_from = Layer),
      by = c("stratum_from", "OmicLayer_from")
    )

  # Join end (to) positions — no need to pull ymin/ymax again
  flow_data <- flow_data |>
    left_join(
      strata_positions |>
        select(group, Layer, end_x = xmin, end_ymin = ymin, end_ymax = ymax) |>
        rename(stratum_to = group, OmicLayer_to = Layer),
      by = c("stratum_to", "OmicLayer_to")
    )

  # Start positions: stack flows
  flow_data <- flow_data |>
    group_by(OmicLayer_from, stratum_from) |> 
    group_split() |> 
    purrr::map(~ {
      .x <- .x |>
        mutate(startdrug_n_fac = as.integer(as.factor(Drug)) - 1)

      total_height <- .x$start_ymax[1] - .x$start_ymin[1]
      total_relative_size <- sum(.x$relative_size)

      .x <- .x |>
        arrange(startdrug_n_fac) |>
        mutate(
          scaled_height = relative_size / total_relative_size * total_height,
          startdrug_ymin = .x$start_ymin[1] + c(0, head(cumsum(scaled_height), -1)),
          startdrug_ymax = startdrug_ymin + scaled_height
        )
      return(.x)
    }) |> 
    bind_rows() |>
    mutate(start_ymin = startdrug_ymin, start_ymax = startdrug_ymax)

  # End positions
flow_data <- flow_data |>
  group_by(OmicLayer_to, stratum_to, stratum_from) |> # Group by OmicLayer_to, stratum_to, AND stratum_from
  group_split() |>
  purrr::map(~ {
    .x <- .x |>
      mutate(enddrug_n_fac = as.integer(as.factor(Drug)) - 1) # Factor for Drug within this stratum_from block

    # Calculate stacking within each stratum_from block
    block_total_height <- sum(.x$relative_size) # Total relative size for this stratum_from block
    
    .x <- .x |>
      arrange(enddrug_n_fac) |> # Arrange by drug within the stratum_from block
      mutate(
        scaled_height_within_stratum = relative_size / block_total_height, # Relative height within its stratum_from block
        end_block_ymin = c(0, head(cumsum(scaled_height_within_stratum), -1)), # Ymin relative to start of its stratum_from block
        end_block_ymax = end_block_ymin + scaled_height_within_stratum # Ymax relative to start of its stratum_from block
      )
    return(.x)
  }) |>
  bind_rows()

# Now, stack the stratum_from blocks within OmicLayer_to, stratum_to
flow_data <- flow_data |>
  group_by(OmicLayer_to, stratum_to) |>
  group_split() |>
  purrr::map(~ {
    # Get unique stratum_from identifiers and their total relative sizes
    stratum_from_summary <- .x |>
      group_by(stratum_from) |>
      summarise(
        stratum_from_relative_size_sum = sum(relative_size),
        # Assuming all flows to this end stratum have the same overall end_ymin/ymax
        stratum_from_overall_ymin = first(end_ymin),
        stratum_from_overall_ymax = first(end_ymax)
      ) |>
      ungroup() |>
      mutate(stratum_from_n_fac = as.integer(as.factor(stratum_from)) - 1) |> # Factor for stratum_from ordering
      arrange(stratum_from_n_fac) # Arrange stratum_from blocks

    total_height_of_dest_stratum <- .x$end_ymax[1] - .x$end_ymin[1]
    total_relative_size_of_dest_stratum <- sum(stratum_from_summary$stratum_from_relative_size_sum)

    stratum_from_summary <- stratum_from_summary |>
      mutate(
        scaled_height_of_stratum_block = stratum_from_relative_size_sum / total_relative_size_of_dest_stratum * total_height_of_dest_stratum,
        stratum_from_base_ymin = .x$end_ymin[1] + c(0, head(cumsum(scaled_height_of_stratum_block), -1))
      )

    # Join the calculated base y-min for each stratum_from block back to the main data
    .x <- .x |>
      left_join(
        stratum_from_summary |> select(stratum_from, stratum_from_base_ymin, scaled_height_of_stratum_block),
        by = "stratum_from"
      ) |>
      mutate(
        # Calculate final end_ymin and end_ymax by adding the block's base Y to its relative Y within the block
        end_ymin = stratum_from_base_ymin + end_block_ymin * scaled_height_of_stratum_block,
        end_ymax = stratum_from_base_ymin + end_block_ymax * scaled_height_of_stratum_block
      )
    return(.x)
  }) |>
  bind_rows()

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
plot_alluvial_from_data <- function(input_data, omics_order, strata_order = NULL, title = NULL) {
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer) # Still loaded, but not directly used for drug_colors here
  library(tidyr)
  library(viridis) # Load the viridis library for scale_fill_viridis_d

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

  strata_positions <- layout_strata_positions(all_strata, strata_order)

  flow_data <- layout_flows_within_strata(input_data, strata_positions)

  # Generate ribbon geometry
  ribbon_data <- generate_alluvium_polygons(flow_data, strata_positions, curve_range = 0.0005)

  # Ensure the 'fill' variable is a factor for discrete palettes like viridis_d
  ribbon_data$fill <- as.factor(ribbon_data$fill)

  strata_positions <- strata_positions %>% mutate(xmid = (xmax + xmin) / 2)
  omics_breaks <- unique(strata_positions$xmid)

  # Build the ggplot object
  final_plot <- ggplot(ribbon_data, aes(x = x, y = y, group = group, fill = fill)) +
    geom_polygon(color = NA, alpha = 0.6) +
    geom_rect(
      data = strata_positions,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), # CORRECTED THIS LINE!
      fill = "white", color = "black", linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = strata_positions,
      aes(x = Layer_Pos, y = (ymin + ymax) / 2, label = group),
      size = 3, hjust = 0.5,
      inherit.aes = FALSE
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5), # Centered title
      axis.title.x = element_text(margin = margin(t = 10), face = "bold"), # Show x-axis title
      axis.text.x = element_text(), # Show x-axis labels
      axis.ticks.x = element_blank(), # Remove x-axis ticks
      plot.margin = margin(t = 10, r = 30, b = 20, l = 10, unit = "pt") # Top, Right, Bottom, Left margins
    ) +
    labs(
      title = title,
      x = "Omics Layer"
    ) +
    scale_x_continuous(breaks = omics_breaks, labels = omics_order) +
    scale_fill_viridis_d(name = "Drug", option = "D", direction = 1) +
    theme(legend.position = "right") # Remove or set to "right" for debugging

  return(final_plot)
}
