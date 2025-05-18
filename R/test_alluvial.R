# Load necessary libraries
install.packages(c('testthat', 'ggplot2', 'dplyr'))
library(testthat)
library(ggplot2)
library(dplyr)
library(vdiffr)
source('./R/utils.R')

# ========== LOAD TEST DATA AND PROCESS FOR TESTING ==========

test_data_path <- 'test-data/multiomics_test_data.csv'
if (!file.exists(test_data_path)) {
  stop(paste("Test data file not found:", test_data_path))
}

mock_input_data <- read.csv(test_data_path)

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

# Layout flows within strata to get the coordinates for alluvia
mock_flow_data_with_coords <- layout_flows_within_strata(mock_input_data, mock_strata_positions)

# Generate alluvium polygons (this is part of the plotting process)
mock_ribbon_data <- generate_alluvium_polygons(mock_flow_data_with_coords, curve_range = 5)

# ========== HELPER FUNCTION FOR TESTING LAYOUT ==========
test_layout_strata <- function(input_data, strata_order) {
  strata_data <- bind_rows(
    input_data %>% select(Layer = OmicLayer_from, group = stratum_from) %>% distinct(),
    input_data %>% select(Layer = OmicLayer_to, group = stratum_to) %>% distinct()
  ) %>%
    mutate(Layer = factor(Layer, levels = correct_omics_order, ordered = TRUE)) %>% # Ensure Layer is ordered factor
    filter(Layer %in% names(strata_order) & group %in% strata_order[[Layer]]) # Filter to relevant groups

  layout_strata_positions(strata_data, strata_order = strata_order)
}

# ========== TESTS ==========

test_that("layout_strata_positions creates position columns", {
  all_strata_test <- bind_rows(
    mock_input_data %>% select(Layer = OmicLayer_from, group = stratum_from) %>% distinct(),
    mock_input_data %>% select(Layer = OmicLayer_to, group = stratum_to) %>% distinct()
  ) %>%
    mutate(Layer = factor(Layer, levels = correct_omics_order, ordered = TRUE)) %>%
    distinct(Layer, group) # Ensure uniqueness here as well
  result <- layout_strata_positions(all_strata_test, strata_order = sample_strata_order)
  expect_true(all(c("xmin", "xmax", "ymin", "ymax", "Layer_Pos") %in% names(result)))
  expect_equal(nrow(result), nrow(unique(all_strata_test[, c("group", "Layer")])))
  expect_false(any(is.na(result$xmin)))
})

test_that("layout_flows_within_strata returns a data frame with flow coordinates", {
  all_strata_test <- bind_rows(
    mock_input_data %>% select(Layer = OmicLayer_from, group = stratum_from) %>% distinct(),
    mock_input_data %>% select(Layer = OmicLayer_to, group = stratum_to) %>% distinct()
  ) %>%
    mutate(Layer = factor(Layer, levels = correct_omics_order, ordered = TRUE)) %>%
    distinct(Layer, group) # Ensure uniqueness here
  mock_strata_positions_test <- layout_strata_positions(all_strata_test, strata_order = sample_strata_order)
  result <- layout_flows_within_strata(mock_input_data, mock_strata_positions_test)
  expect_true(all(c("start_x", "end_x", "start_ymin", "end_ymin", "start_ymax", "end_ymax") %in% names(result)))
  expect_equal(nrow(result), nrow(mock_input_data))
})

test_that("generate_alluvium_polygons returns valid polygon data", {
  all_strata_test <- bind_rows(
    mock_input_data %>% select(Layer = OmicLayer_from, group = stratum_from) %>% distinct(),
    mock_input_data %>% select(Layer = OmicLayer_to, group = stratum_to) %>% distinct()
  ) %>%
    mutate(Layer = factor(Layer, levels = correct_omics_order, ordered = TRUE)) %>%
    distinct(Layer, group) # Ensure uniqueness here
  mock_strata_positions_test <- layout_strata_positions(all_strata_test, strata_order = sample_strata_order)
  mock_flow_data_with_coords_test <- layout_flows_within_strata(mock_input_data, mock_strata_positions_test)
  mock_ribbon_data_test <- generate_alluvium_polygons(mock_flow_data_with_coords_test, n=50, curve_range = 5)
  expect_true(all(c("x", "y", "group", "fill") %in% names(mock_ribbon_data_test)))
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

test_that("plot_alluvial_from_data uses the correct number of ribbon groups", {
  plot_output <- plot_alluvial_from_data(mock_input_data, omics_order = correct_omics_order, strata_order = sample_strata_order)
  ribbon_data <- ggplot_build(plot_output)$data[[1]]
  expected_number_of_groups <- nrow(mock_input_data)
  expect_equal(length(unique(ribbon_data$group)), expected_number_of_groups)
})

test_that("plot_alluvial_from_data uses a manual color scale for drugs", {
  plot_output <- plot_alluvial_from_data(mock_input_data, omics_order = correct_omics_order, strata_order = sample_strata_order)
  scales_list <- plot_output$scales$scales
  fill_scale <- NULL
  for (scale in scales_list) {
    if ("fill" %in% scale$aesthetics) { # Check if "fill" is one of the aesthetics
      fill_scale <- scale
      break
    }
  }
  expect_false(is.null(fill_scale), info = "Fill scale not found")
  if (!is.null(fill_scale)) {
    expect_s3_class(fill_scale, "ScaleDiscrete")
    scale_breaks <- fill_scale$get_breaks()
    unique_drugs <- unique(mock_input_data$Drug)
    expect_equal(length(scale_breaks), length(unique_drugs))
  }
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
omics_order_for_plot <- c("Genomics", "Transcriptomics", "Proteomics", "Metabolomics", "Lipidomics")

try(print(plot_alluvial_from_data(mock_input_data, omics_order = omics_order_for_plot)))
