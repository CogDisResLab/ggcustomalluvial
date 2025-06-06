# Load necessary libraries
library(testthat)
library(ggplot2)
library(tidyverse)
library(purrr)
library(dplyr)
source('./R/utils.R')

# ========== LOAD TEST DATA AND PROCESS FOR TESTING ==========

test_data_path <- 'test-data/multiomics_test_data.csv'
if (!file.exists(test_data_path)) {
  stop(paste("Test data file not found:", test_data_path))
}

mock_input_data <- read.csv(test_data_path)

actual_input_data <- read.csv('test-data/multiomics_test_data.csv')

# Define a correct omics order
correct_omics_order <- c("Genomics", "Transcriptomics", "Proteomics", "Metabolomics", "Lipidomics")

# Ensure the 'OmicLayer_from' and 'OmicLayer_to' columns are factors with the correct order
mock_input_data <- mock_input_data %>%
  mutate(
    OmicLayer_from = factor(OmicLayer_from, levels = correct_omics_order, ordered = TRUE),
    OmicLayer_to = factor(OmicLayer_to, levels = correct_omics_order, ordered = TRUE)
  )

# Collect all unique strata from the input data
all_strata <- bind_rows(
  mock_input_data %>% select(Layer = OmicLayer_from, group = stratum_from) %>% distinct(),
  mock_input_data %>% select(Layer = OmicLayer_to, group = stratum_to) %>% distinct()
) %>% distinct()

# Ensure the 'Layer' column is a factor with the correct order
all_strata <- all_strata %>%
  mutate(Layer = factor(Layer, levels = correct_omics_order, ordered = TRUE))

# Define a sample strata order for testing
sample_strata_order <- list(
  "Genomics" = unique(all_strata$group[all_strata$Layer == "Genomics"]),
  "Transcriptomics" = unique(all_strata$group[all_strata$Layer == "Transcriptomics"]),
  "Proteomics" = unique(all_strata$group[all_strata$Layer == "Proteomics"]),
  "Metabolomics" = unique(all_strata$group[all_strata$Layer == "Metabolomics"]),
  "Lipidomics" = unique(all_strata$group[all_strata$Layer == "Lipidomics"])
)

# Layout strata positions based on the unique strata
mock_strata_positions <- layout_strata_positions(all_strata, strata_order = sample_strata_order)

# ========== HELPER FUNCTION FOR TESTING LAYOUT ==========
test_layout_strata <- function(input_data, strata_order) {
  strata_data <- bind_rows(
    input_data %>% select(Layer = OmicLayer_from, group = stratum_from) %>% distinct(),
    input_data %>% select(Layer = OmicLayer_to, group = stratum_to) %>% distinct()
  ) %>%
    mutate(Layer = factor(Layer, levels = correct_omics_order, ordered = TRUE)) %>% # Ensure Layer is ordered factor
    filter(Layer %in% names(strata_order) & group %in% strata_order[[Layer]]) # Filter to relevant groups

  return(layout_strata_positions(strata_data, strata_order = strata_order))
}

test_position_adjacency <- function(positions, yvar_1, yvar_2) {
  truth_table <- c()
  positions <- round(positions[, c(yvar_1, yvar_2)], 6)
  for (i in 2:nrow(positions)) {
    prev_vector <- positions[i - 1, ]
    vector <- positions[i, ]
    if (any(is.na(prev_vector)) || any(is.na(vector))) {
      truth_table <- c(truth_table, FALSE)
    } else {
      prev_vector <- as.vector(unlist(prev_vector))
      vector <- as.vector(unlist(vector))
      truth_table <- c(truth_table, prev_vector[2] == vector[1])
   }
  }
  return(truth_table)
}

test_split_positions_by_size <- function(input_data) {
  truth_table <- c()

  # Group by OmicLayer_from, stratum_from, and Drug
  grouped_input <- input_data %>%
    group_by(OmicLayer_from, stratum_from, Drug) %>%
    group_split()

  for (group in grouped_input) {
    positions <- split_positions_by_size(group)
    if (nrow(positions) < 2) {
      next  # Skip groups with only one entry
    }

  truth_table <- c(truth_table, test_position_adjacency(positions))
  }

  return(truth_table)
}

generate_input_for_flows_within_function <- function(input_data, omics_order=correct_omics_order, strata_order=NULL) {
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
  return(input_data)
}

generate_input_for_polygons_function <- function(input_data, omics_order=correct_omics_order, strata_order=NULL) {
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
  if (!("Drug" %in% colnames(flow_data))) {
    stop("Drug not in colnames of flow_data from layout_flows_within_strata")
  }
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
  return(flow_data)
}

generate_input_for_flow_function <- function(input_data, omics_order=correct_omics_order, strata_order=NULL) {
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

  flow_data <- input_data
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
        mutate(total_size = sum(size)) %>%
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
        )
      processed_list_1[[length(processed_list_1) + 1]] <- df_processed
    }
  }
  # Combine the processed sub-dataframes from list_1 back into one
  processed_flow_data_1 <- bind_rows(processed_list_1)
  return(processed_flow_data_1)
}

# ========== TESTS ==========
test_that("split_positions_by_size returns adjacent positions", {
  split_input_data <- generate_input_for_flow_function(mock_input_data)
  expect_true(all(test_split_positions_by_size(split_input_data)))
})

test_that("layout_strata_positions creates position columns", {
  all_strata_test <- bind_rows(
    mock_input_data %>% select(Layer = OmicLayer_from, group = stratum_from) %>% distinct(),
    mock_input_data %>% select(Layer = OmicLayer_to, group = stratum_to) %>% distinct()
  ) %>%
    mutate(Layer = factor(Layer, levels = correct_omics_order, ordered = TRUE)) %>%
    distinct(Layer, group) # Ensure uniqueness here as well
  result <- layout_strata_positions(all_strata_test, strata_order = sample_strata_order)
  result_adjacent <- result |> group_by(Layer) |> as.data.frame()
  result_adjacent <- result |> group_split(Layer) |> purrr::map(function(X) test_position_adjacency(X, 'ymin', 'ymax')) |> map(function(X) expect_true(all(X)))
  expect_true(all(c("xmin", "xmax", "ymin", "ymax", "Layer_Pos") %in% names(result)))
  expect_equal(nrow(result), nrow(unique(all_strata_test[, c("group", "Layer")])))
  expect_false(any(is.na(result$xmin)))
})

test_that("layout_flows_within_strata returns a data frame with flow coordinates", {
  split_input_data <- generate_input_for_flows_within_function(mock_input_data)
  result <- layout_flows_within_strata(split_input_data, mock_strata_positions)
  n_result <- layout_flows_within_strata(split_input_data, mock_strata_positions)
  expect_true(all(c("start_x", "end_x", "start_ymin", "end_ymin", "start_ymax", "end_ymax") %in% names(result)))
  expect_equal(nrow(result), nrow(mock_input_data))
})

test_that("generate_alluvium_polygons returns valid polygon data", {
  flow_data <- layout_flows_within_strata(mock_input_data, mock_strata_positions)
  mock_ribbon_data_test <- generate_alluvium_polygons(flow_data, mock_strata_positions, n=50, curve_range = 5)
  expect_gt(nrow(mock_ribbon_data_test), 0)
})

test_that("plot_alluvial_from_data returns a ggplot object", {
  plot_output <- plot_alluvial_from_data(mock_input_data, omics_order = correct_omics_order, strata_order = sample_strata_order)
  expect_s3_class(plot_output, "ggplot")
})

test_that("plot_alluvial_from_data throws error if omics_order is NULL", {
  expect_error(plot_alluvial_from_data(mock_input_data, omics_order = NULL),
               regexp = "Error: The 'omics_order' argument must be provided.")
})

test_that("plot_alluvial_from_data with correct omics_order does not throw error", {
  expect_silent(plot_alluvial_from_data(mock_input_data, omics_order = correct_omics_order, strata_order = sample_strata_order))
})

test_that("plot_alluvial_from_data throws error for incorrect flow direction with omics_order", {
  incorrect_omics_order_reversed <- rev(correct_omics_order)
  expect_error(plot_alluvial_from_data(mock_input_data, omics_order = incorrect_omics_order_reversed),
               regexp = "Error: Some flows are moving against or staying within the specified omics layer order.")

  incorrect_omics_order_mixed <- c("Genomics", "Proteomics", "Transcriptomics", "Metabolomics", "Lipidomics")
  expect_error(plot_alluvial_from_data(mock_input_data, omics_order = incorrect_omics_order_mixed),
               regexp = "Error: Some flows are moving against or staying within the specified omics layer order.")
})

test_that("plot_alluvial_from_data generates ribbon data with correct structure", {
  plot_output <- plot_alluvial_from_data(mock_input_data, omics_order = correct_omics_order, strata_order = sample_strata_order)
  ribbon_data <- ggplot_build(plot_output)$data[[1]]
  expect_true(all(c("x", "y", "group", "fill") %in% names(ribbon_data)))
  expect_gt(nrow(ribbon_data), 0)
})

test_that("layout_strata_positions correctly handles layers with few strata", {
  # Create a subset of the data where Lipidomics has only 2 strata
  few_strata_data <- mock_input_data %>%
    filter(OmicLayer_to == "Lipidomics" & stratum_to %in% c("RegionA", "RegionB")) %>%
    bind_rows(mock_input_data %>% filter(OmicLayer_to != "Lipidomics")) %>%
    mutate(
      OmicLayer_from = factor(OmicLayer_from, levels = correct_omics_order, ordered = TRUE),
      OmicLayer_to = factor(OmicLayer_to, levels = correct_omics_order, ordered = TRUE)
    )

  # Define a strata order for this subset
  few_strata_order <- list(
    "Genomics" = unique(few_strata_data$group[few_strata_data$Layer == "Genomics"]),
    "Transcriptomics" = unique(few_strata_data$group[few_strata_data$Layer == "Transcriptomics"]),
    "Proteomics" = unique(few_strata_data$group[few_strata_data$Layer == "Proteomics"]),
    "Metabolomics" = unique(few_strata_data$group[few_strata_data$Layer == "Metabolomics"]),
    "Lipidomics" = c("RegionB", "RegionA") # Specific order for Lipidomics
  )

  # Create strata_data locally within the test, selecting only Layer and group IMMEDIATELY after bind_rows
  strata_data_few_original <- bind_rows(
    few_strata_data %>% select(Layer = OmicLayer_from, group = stratum_from) %>% distinct(),
    few_strata_data %>% select(Layer = OmicLayer_to, group = stratum_to) %>% distinct()
  ) %>%
    mutate(Layer = factor(Layer, levels = correct_omics_order, ordered = TRUE)) %>%
    distinct(Layer, group) %>%
    select(Layer, group) # Ensure only Layer and group are present RIGHT HERE

  # Pass a copy to layout_strata_positions
  layout_result <- layout_strata_positions(strata_data_few_original %>% tibble::as_tibble() %>% as.data.frame(),
                                          strata_order = few_strata_order)

  # Find the ymin and ymax for Lipidomics strata
  lipidomics_strata <- layout_result %>%
    filter(Layer == "Lipidomics") %>%
    arrange(ymin) # Order by ymin to match the top-to-bottom layout

  ymin_lipidomics <- lipidomics_strata$ymin
  ymax_lipidomics <- lipidomics_strata$ymax

  # Check if the y-range for Lipidomics spans approximately 0 to 1
  expect_true(min(ymin_lipidomics) < 0.1) # Should be close to 0
  expect_true(max(ymax_lipidomics) > 0.9) # Should be close to 1

  # Add a specific test for strata order in Lipidomics
  expect_equal(lipidomics_strata$group[1], "RegionB")
  expect_equal(lipidomics_strata$group[2], "RegionA")
})

print(plot_alluvial_from_data(mock_input_data, omics_order = correct_omics_order, strata_order = sample_strata_order))
