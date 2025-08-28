# ggcustomalluvial: A ggplot2 Extension for Alluvial Diagrams
`ggcustomalluvial` is a powerful R package that extends the ggplot2 framework to create beautiful and informative alluvial diagrams. This package is specifically designed to handle complex data, such as multi-omics datasets, and allows for custom ordering and styling of flows and strata.

The package handles all the necessary data manipulation and layout calculations, freeing you to focus on visualizing the relationships and flows within your data.

## Features
Integrated with ggplot2: Use ggcustomalluvial seamlessly with the ggplot2 ecosystem.

Custom Layout: Define the precise order of your data layers and strata to create meaningful visualizations.

Multi-omics Ready: Optimized to handle multi-layered data, such as genomics, transcriptomics, and proteomics.

Clean & Descriptive Output: The package produces publication-quality plots that are easy to interpret.

## Installation
You can install the package directly from GitHub using the devtools package:
```r
install.packages("devtools")
devtools::install_github("CogDisResLab/ggcustomalluvial")
```


## Usage
The core function of the package is `plot_alluvial_from_data()`. You provide your data frame, a specified order for your omics layers, and optionally, a custom order for the strata within each layer.

```r
# Load the package and example data
library(ggcustomalluvial)
library(dplyr)
library(ggplot2)

# Assuming 'multiomics_test_data.csv' is available
# Here's a quick example of a sample data structure
# that the function expects.

df <- data.frame(
  OmicLayer_from = c("Genomics", "Genomics", "Transcriptomics"),
  stratum_from = c("Tumor A", "Tumor B", "Tumor A"),
  OmicLayer_to = c("Transcriptomics", "Transcriptomics", "Proteomics"),
  stratum_to = c("Tumor A", "Tumor B", "Tumor A"),
  value = c(10, 5, 8),
  Drug = c("DrugX", "DrugY", "DrugX")
)

# Define the order of omics layers
omics_order <- c("Genomics", "Transcriptomics", "Proteomics", "Metabolomics", "Lipidomics")

# You can also define a specific order for strata within each layer
strata_order <- list(
  Genomics = c("Tumor A", "Tumor B"),
  Transcriptomics = c("Tumor A", "Tumor B"),
  Proteomics = c("Tumor A"),
  Metabolomics = c(),
  Lipidomics = c()
)

# Generate and plot the alluvial diagram
alluvial_plot <- plot_alluvial_from_data(
  data = df,
  omics_order = omics_order,
  strata_order = strata_order,
  fill_by = "Drug" # Optional: fill by a specific column
)

print(alluvial_plot)

# To save the plot
# ggsave("alluvial_plot.png", plot = alluvial_plot)
```
Contributing
We welcome contributions! Please feel free to submit a pull request or open an issue if you encounter any bugs or have suggestions for new features.

License
This project is licensed under the MIT License.
