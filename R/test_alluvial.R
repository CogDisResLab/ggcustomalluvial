library(ggplot2)
library(dplyr)

# Load data
df_alluvial <- read.csv('test_data/alluvial_long.csv')

# Plot using the all-in-one wrapper
final_plot <- plot_alluvial_from_data(df_alluvial)

# Display
print(final_plot)
