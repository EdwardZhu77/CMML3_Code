# --- Step2_Monitoring_Equilibration.R ---
# Purpose: Generate equilibration plots (Density, Potential Energy, Pressure)
#          for a single, specified simulation run.
# WARNING: This script assumes input files exist and are correctly formatted.

# --- 0. Load necessary packages ---
library(tidyverse)
library(fs)

# --- 1. User Settings ---

# Base directory where MD_Result is located (assuming script is run from project base)
base_dir <- getwd()
# Specific simulation run directory to analyze
simulation_folder <- path(base_dir, "MD_Result", "280K", "280K_2")
# Label for the plot title
temperature_label <- "280K (Rep 2)"

# Define metrics to plot and their properties
metrics_for_single_plot <- list(
  `Density` = list(
    filename = "density.xvg",
    ylabel = "Density (kg/m³)",
    needs_conversion = TRUE # Assuming density.xvg is from NPT, time in ps
  ),
  `Potential Energy` = list(
    filename = "potential.xvg",
    ylabel = "Potential Energy (kJ/mol)",
    needs_conversion = TRUE # Assuming potential.xvg is from EM, time in ps
  ),
  `Pressure` = list(
    filename = "pressure.xvg",
    ylabel = "Pressure (bar)",
    needs_conversion = TRUE # Assuming pressure.xvg is from NPT, time in ps
  )
)

# Plotting aesthetics
plot_color <- "#6495ED" # CornflowerBlue for 280K
time_ps_to_ns_factor <- 0.001 # Conversion factor

# Output directory and filename
output_dir <- path(simulation_folder) # Save plot in the simulation folder
safe_temp_label <- gsub("[^a-zA-Z0-9_]", "_", temperature_label)
output_filename <- path(output_dir, paste0("equilibration_monitor_", safe_temp_label, ".png"))

# Time range for calculating average values (adjust if needed)
avg_start_time_ns <- 0.1 # Start averaging after 0.1 ns

# --- 2. Read Data ---
all_single_sim_data <- list() # Initialize list to store data frames

# Loop through each metric defined above
for (metric_label in names(metrics_for_single_plot)) {
  metric_details <- metrics_for_single_plot[[metric_label]]
  file_path <- path(simulation_folder, metric_details$filename)
  message("  Reading: ", file_path)

  # Read lines, skipping comments/headers
  lines <- readLines(file_path)
  data_lines <- lines[!startsWith(lines, "@") & !startsWith(lines, "#") & nzchar(trimws(lines))]

  # Convert data lines to a data frame
  temp_df <- read.table(textConnection(data_lines), header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  current_data_df <- temp_df[, 1:2] # Select first two columns
  colnames(current_data_df) <- c("time_raw", "value") # Assign column names

  # Convert time units if needed and add metric label
  current_data_df <- current_data_df %>%
    filter(!is.na(time_raw) & !is.na(value)) %>% # Remove any NA rows
    mutate(
      time_ns = if (metric_details$needs_conversion) time_raw * time_ps_to_ns_factor else time_raw, # Convert ps to ns
      metric_type = factor(metric_label, levels = names(metrics_for_single_plot)) # Add metric type factor
    ) %>%
    select(time_ns, value, metric_type) # Keep relevant columns

  # Add the processed data frame to the list
  all_single_sim_data[[metric_label]] <- current_data_df
  message("  Processed ", nrow(current_data_df), " data points for ", metric_label)
}

# --- 3. Combine Data ---
# Combine data frames for all metrics into one
combined_single_sim_data <- bind_rows(all_single_sim_data)

# Get the maximum time for plotting limits and text placement
max_time_ns <- max(combined_single_sim_data$time_ns, na.rm = TRUE)

# --- 4. Calculate Averages ---
# Calculate average density and pressure after the defined start time
message("\n--- Calculating average Density and Pressure (", avg_start_time_ns, "ns - ", max_time_ns, "ns) ---")
average_values <- combined_single_sim_data %>%
  filter(metric_type %in% c("Density", "Pressure")) %>% # Select only Density and Pressure
  filter(time_ns >= avg_start_time_ns & time_ns <= max_time_ns) %>% # Filter by time range
  group_by(metric_type) %>% # Group by metric
  summarise(avg_value = mean(value, na.rm = TRUE), .groups = 'drop') %>% # Calculate mean
  mutate(text_x_pos = max_time_ns * 0.85) # Define X position for text label

print("Calculated average values:")
print(average_values)

# --- 5. Create Combined Plot ---
# Create labels for facets based on ylabel defined earlier
facet_labels <- sapply(names(metrics_for_single_plot), function(name) {
  metrics_for_single_plot[[name]]$ylabel
})

# Generate the plot using ggplot2
combined_plot <- ggplot(data = combined_single_sim_data,
                        mapping = aes(x = time_ns, y = value)) +
  # Plot the time series lines
  geom_line(color = plot_color, alpha = 0.9, linewidth = 0.6) +
  # Add horizontal dashed lines for averages
  geom_hline(data = average_values,
             aes(yintercept = avg_value),
             linetype = "dashed", color = "darkblue", linewidth = 0.8) +
  # Add text labels for the average values
  geom_text(data = average_values,
            aes(label = sprintf("Avg: %.2f", avg_value), y = avg_value, x = text_x_pos),
            hjust = 1, vjust = -1.0, color = "darkblue", size = 3.5, fontface = "bold") +
  # Create separate panels (facets) for each metric
  facet_wrap(~ metric_type, ncol = 1, scales = "free_y", strip.position = "left",
             labeller = labeller(metric_type = facet_labels)) +
  # Set plot titles and axis labels
  labs(
    title = paste("Basic Equilibration Metrics for", temperature_label),
    x = "Time (ns)",
    y = NULL # Y axis label handled by facet strips
  ) +
  # Apply theme settings
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0.5, "lines")
  )

# --- 6. Save Plot ---
# Save the generated plot to a PNG file
ggsave(
  filename = output_filename,
  plot = combined_plot,
  width = 8, height = 9, dpi = 300
)
message("Combined equilibration plot saved to: ", output_filename)
