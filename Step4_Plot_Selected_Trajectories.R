# --- Step4_Plot_Selected_Trajectories_Manual.R ---
# Purpose: Plot RMSD (vs initial, vs crystal) and Rg for a MANUALLY
#          specified replicate for each temperature.
# WARNING: Assumes specified replicate XVG files exist.

# --- 0. Load necessary packages ---
library(tidyverse)
library(fs)
# library(RColorBrewer) # Optional for colors

# --- 1. User Settings ---

# Base directory where MD_Result and Analysis_Output are located
base_dir <- getwd()
results_base_dir <- path(base_dir, "MD_Result")
analysis_output_dir <- path(base_dir, "Analysis_Output") # Where output plot goes
dir_create(analysis_output_dir) # Create output directory if needed

# Temperatures to plot
temperatures <- c("280K", "300K", "320K")

# --- !!! MANUAL SPECIFICATION OF REPLICATES !!! ---
# Edit the numbers in `replicate_to_plot` to choose which replicate (1, 2, or 3)
# to plot for each corresponding temperature in the `temperatures` vector.
# Example: c(2, 3, 1) means use Rep 2 for 280K, Rep 3 for 300K, Rep 1 for 320K.
selected_replicate_numbers <- c(2, 3, 1)
# --- End Manual Specification ---

# Create a tibble mapping temperature to the chosen replicate number
if(length(temperatures) != length(selected_replicate_numbers)) {
  stop("Error: Length of 'temperatures' and 'selected_replicate_numbers' must match.")
}
manual_selection_info <- tibble(
  temperature = factor(temperatures, levels = temperatures),
  replicate_id_num = selected_replicate_numbers
)

print("--- Using MANUALLY specified replicates for plotting: ---")
print(manual_selection_info)

# --- 1.1 Define metrics to plot ---
metrics_to_plot <- list(
  `RMSD (Initial Ref)` = list( # RMSD vs Minimized Structure (em.tpr)
    key = "rmsd_xtal",
    filename_base = "rmsd_xtal", # Base name without _run_name.xvg
    needs_conversion = FALSE # Assuming gmx rms -tu ns
  ),
  `RMSD (MD Start Ref)` = list( # RMSD vs Production Start Structure (md*.tpr)
    key = "rmsd",
    filename_base = "rmsd",      # Base name without _run_name.xvg
    needs_conversion = FALSE # Assuming gmx rms -tu ns
  ),
  `Rg` = list(
    key = "gyrate",
    filename_base = "gyrate",    # Base name without _run_name.xvg
    needs_conversion = TRUE  # Assuming gmx gyrate output time is ps
  )
)

# --- 1.2 Color scheme ---
temp_color_map <- c(
  "280K" = "#6495ED",  # CornflowerBlue
  "300K" = "#3CB371",  # MediumSeaGreen
  "320K" = "#CD5C5C"   # IndianRed
)

# --- 1.3 Other settings ---
time_ps_to_ns_factor <- 0.001
output_plot_filename <- path(analysis_output_dir, "Manually_Selected_Trajectories_Plot.png") # Changed output name

# --- 2. Read Data for Manually Selected Replicates ---
message("\n--- Reading RMSD and Rg data ONLY for MANUALLY selected replicates ---")
all_selected_data_list <- list() # Initialize list

# Loop through each METRIC to plot
for (metric_label in names(metrics_to_plot)) {
  metric_details <- metrics_to_plot[[metric_label]]
  message("Processing Metric: ", metric_label)

  # Loop through each ROW in the manual_selection_info table
  for (row_idx in 1:nrow(manual_selection_info)) {
    temp <- manual_selection_info$temperature[row_idx]
    rep_id_num <- manual_selection_info$replicate_id_num[row_idx] # Get the chosen replicate number

    # Construct the run ID (e.g., "280K_2")
    run_id <- paste0(temp, "_", rep_id_num)

    # Construct the specific XVG filename and path
    filename <- paste0(metric_details$filename_base, "_", run_id, ".xvg")
    file_path <- path(results_base_dir, temp, run_id, filename)

    message("  Reading: ", file_path)
    # Read lines, skip headers/comments
    lines <- readLines(file_path)
    data_lines <- lines[!startsWith(lines, "@") & !startsWith(lines, "#") & nzchar(trimws(lines))]

    # Process into data frame
    temp_df <- read.table(textConnection(data_lines), header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
    current_data_df <- temp_df[, 1:2]
    colnames(current_data_df) <- c("time_raw", "value")

    # Mutate and select relevant columns
    current_data_df <- current_data_df %>%
      filter(!is.na(time_raw) & !is.na(value)) %>%
      mutate(
        time_ns = if (metric_details$needs_conversion) time_raw * time_ps_to_ns_factor else time_raw,
        temperature = factor(temp, levels = temperatures), # Ensure temperature is a factor
        metric_type = factor(metric_label, levels = names(metrics_to_plot)) # Add metric type factor
      ) %>%
      select(time_ns, value, temperature, metric_type)

    # Add the processed data frame to the list
    all_selected_data_list[[paste(metric_label, run_id, sep = "_")]] <- current_data_df
    message("  Processed ", nrow(current_data_df), " data points for ", metric_label, " from ", run_id)

  } # End loop through selected runs
} # End loop through metrics

# Combine all data frames into one
combined_selected_data <- bind_rows(all_selected_data_list)

# --- 3. Generate Faceted Plot ---
message("\n--- Generating faceted plot for manually selected trajectories ---")

# Generate the plot
selected_traj_plot <- ggplot(data = combined_selected_data,
                             mapping = aes(x = time_ns, y = value, color = temperature)) +
  geom_line(alpha = 0.9, linewidth = 0.7) +
  scale_color_manual(values = temp_color_map, name = "Temperature") +
  facet_wrap(~ metric_type, scales = "free_y", ncol = 1) +
  labs(
    title = "RMSD & Rg for Manually Selected Replicates", # Updated title
    x = "Time (ns)",
    y = "Value (units vary by plot)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold", size = 11)
  )

# --- 4. Save Plot ---
ggsave(
  filename = output_plot_filename,
  plot = selected_traj_plot,
  width = 8, height = 10, dpi = 300
)
message("Faceted plot for manually selected trajectories saved to: ", output_plot_filename)

message("\n--- R Script Finished ---")