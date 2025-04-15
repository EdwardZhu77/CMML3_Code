# --- Step3_Plot_Select_Rg.R ---
# Purpose: 1. Plot average Radius of Gyration (Rg) across replicates per temperature.
#          2. Calculate Rg variance for each replicate after equilibration.
#          3. Identify and report the most stable replicate (lowest variance) per temperature.
# WARNING: Assumes input files exist and are correctly formatted.

# --- 0. Load necessary packages ---
library(tidyverse)
library(fs)

# --- 1. User Settings ---

# Base directory where MD_Result is located
base_dir <- getwd()
results_base_dir <- path(base_dir, "MD_Result") # Path to the main results folder

# Simulation parameters
temperatures <- c("280K", "300K", "320K")
num_replicates <- 3
replicate_pattern <- "{temp}_{rep_num}" # Folder name pattern within TEMP folders

# Rg file details
metric_key <- "gyrate"
metric_filename_pattern <- "gyrate_{run_name}.xvg" # Assumes filename includes run_name like gyrate_280K_1.xvg
metric_ylabel <- "Rg (nm)"
metric_title_base <- "Radius of Gyration (Rg) vs Time"

# Plotting aesthetics
color_map <- c(
  "280K" = "#6495ED",  # CornflowerBlue
  "300K" = "#3CB371",  # MediumSeaGreen
  "320K" = "#CD5C5C"   # IndianRed
)

# Time conversion and output directory
time_unit_factor <- 0.001 # ps to ns
output_dir <- path(base_dir, "Analysis_Output") # Save plots/results here
dir_create(output_dir) # Create output directory if it doesn't exist

# Equilibration time for variance calculation
equilibration_ns <- 10 # ns

# --- 2. Read all Rg data ---
all_data_list <- list() # Initialize list

# Loop through temperatures and replicates to find Rg files
for (temp in temperatures) {
  message("Processing Temperature: ", temp)
  for (i in 1:num_replicates) {
    run_name <- gsub("\\{temp\\}", temp, gsub("\\{rep_num\\}", i, replicate_pattern))
    filename <- gsub("\\{run_name\\}", run_name, metric_filename_pattern)
    file_path <- path(results_base_dir, temp, run_name, filename) # Path to the specific Rg file

    message("  Reading: ", file_path)
    # Read lines, skip headers/comments
    lines <- readLines(file_path)
    data_lines <- lines[!startsWith(lines, "@") & !startsWith(lines, "#") & nzchar(trimws(lines))]

    # Process data lines into a data frame
    temp_df <- read.table(textConnection(data_lines), header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
    current_data_df <- temp_df[, 1:2]
    colnames(current_data_df) <- c("time_ps", "value")

    # Mutate and select relevant columns
    current_data_df <- current_data_df %>%
      filter(!is.na(time_ps) & !is.na(value)) %>%
      mutate(
        time_ns = time_ps * time_unit_factor,
        temperature = factor(temp, levels = temperatures),
        replicate_id = factor(i) # Keep replicate ID
      ) %>%
      select(time_ns, value, temperature, replicate_id)

    # Add the processed data frame to the list
    all_data_list[[run_name]] <- current_data_df # Use run_name as list index
    message("  Processed ", nrow(current_data_df), " data points for ", run_name)
  }
}

# Combine all data frames into one
combined_data <- bind_rows(all_data_list, .id = "run_id") # .id adds run_name as a column

# --- 3. Calculate Average Rg and Standard Deviation per Time Point ---
# Group by temperature and time, calculate mean/sd across replicates
average_rg_data <- combined_data %>%
  group_by(temperature, time_ns) %>%
  summarise(
    avg_rg = mean(value, na.rm = TRUE),
    sd_rg = sd(value, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Calculate ribbon limits (handle potential NAs if only 1 rep at a time point)
  mutate(
    sd_rg = ifelse(is.na(sd_rg), 0, sd_rg),
    ymin = avg_rg - sd_rg,
    ymax = avg_rg + sd_rg
  )

# --- 4. Plot Average Rg with Standard Deviation Ribbon ---
plot_avg_rg <- ggplot(data = average_rg_data,
                      mapping = aes(x = time_ns, color = temperature, fill = temperature)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, linetype = "blank") + # SD ribbon
  geom_line(aes(y = avg_rg), linewidth = 0.8) + # Mean line
  scale_color_manual(values = color_map, name = "Temperature") + # Colors for lines
  scale_fill_manual(values = color_map, name = "Std Dev Range") + # Colors for ribbons
  labs(
    title = paste(metric_title_base, "- Average across Replicates"),
    x = "Time (ns)",
    y = metric_ylabel
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# Save the average Rg plot
output_filename_avg <- path(output_dir, "Rg_Average_vs_Time.png")
ggsave(
  filename = output_filename_avg,
  plot = plot_avg_rg,
  width = 10, height = 6, dpi = 300
)
message("Average Rg plot saved to: ", output_filename_avg)

# --- 5. Calculate Fluctuation (Variance) per Replicate ---
# Calculate variance for each individual run after equilibration time
variance_per_replicate <- combined_data %>%
  filter(time_ns >= equilibration_ns) %>% # Consider only production phase
  group_by(temperature, replicate_id, run_id) %>% # Group by each specific run
  summarise(
    rg_variance = var(value, na.rm = TRUE), # Calculate variance
    rg_mean_prod = mean(value, na.rm = TRUE), # Also get mean for reference
    .groups = 'drop' # Ungroup
  )

# --- 6. Identify Most Stable Replicate per Temperature ---
# For each temperature, find the replicate with the minimum variance
most_stable_replicates <- variance_per_replicate %>%
  group_by(temperature) %>% # Group by temperature
  arrange(rg_variance, .by_group = TRUE) %>% # Arrange by variance within each group
  slice_min(order_by = rg_variance, n = 1) %>% # Select the row with minimum variance
  ungroup() # Ungroup

# --- 7. Report Results ---
print("--- Rg Variance Calculation per Replicate (Production Phase) ---")
print(variance_per_replicate)

print("--- Most Stable Replicate per Temperature (Lowest Rg Variance) ---")
print(most_stable_replicates)

# Save the table of most stable replicates
output_table_stable <- path(output_dir, "most_stable_replicates_Rg_variance.csv")
write_csv(most_stable_replicates, output_table_stable)
message("Table of most stable replicates saved to: ", output_table_stable)

message("\n--- R Script Finished ---")