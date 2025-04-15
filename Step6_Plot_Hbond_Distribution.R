# --- Step5b_Plot_Hbond_Distribution.R ---
# Purpose: Read H-bond data calculated by Step5a and generate a
#          box plot with significance annotations. Performs statistical tests.
# WARNING: Assumes required XVG files exist in the Compare_Selected directory.
#          Requires ggpubr, car, dunn.test packages for stats/plotting.

# --- 0. Load necessary packages ---
library(tidyverse)
library(fs)
library(ggpubr)    # For stat_compare_means
# library(car)       # For Levene's test (Optional: can comment out if not needed)
# library(dunn.test) # For Dunn's test (Optional: can comment out if not needed)

# --- 1. User Settings ---

# Base directory where Compare_Selected and Analysis_Output are located
base_dir <- getwd()
compare_dir <- path(base_dir, "Compare_Selected") # Directory with hbond_num_*.xvg files
analysis_output_dir <- path(base_dir, "Analysis_Output") # Where plot will be saved
dir_create(analysis_output_dir) # Ensure output dir exists

# Temperatures to analyze
temperatures <- c("280K", "300K", "320K")
hbond_filename_pattern <- "hbond_num_{temp}.xvg"

# Plotting aesthetics
temp_color_map <- c(
  "280K" = "#6495ED",  # CornflowerBlue
  "300K" = "#3CB371",  # MediumSeaGreen
  "320K" = "#CD5C5C"   # IndianRed
)

# Time conversion (though time column isn't strictly needed for boxplot)
time_ps_to_ns_factor <- 0.001

# Output filename
output_plot_filename <- path(analysis_output_dir, "Hbond_Distribution_Plot.png")

# --- 2. Read and process H-bond data ---
message("\n--- Reading and processing hydrogen bond data from: ", compare_dir, " ---")
all_hbonds_list <- list() # Initialize list

# Loop through temperatures to find files
for (temp in temperatures) {
  filename <- gsub("\\{temp\\}", temp, hbond_filename_pattern)
  file_path <- path(compare_dir, filename)
  message("Processing: ", file_path)

  # Read lines, skip headers/comments
  lines <- readLines(file_path)
  data_lines <- lines[!startsWith(lines, "@") & !startsWith(lines, "#") & nzchar(trimws(lines))]

  # Process into data frame
  temp_df <- read.table(textConnection(data_lines), header = FALSE, col.names = c("time_ps", "hbonds"))

  # Mutate and select columns
  current_hbond_data <- temp_df %>%
    filter(!is.na(time_ps) & !is.na(hbonds)) %>%
    mutate(
      time_ns = time_ps * time_ps_to_ns_factor,
      temperature = factor(temp, levels = temperatures)
    ) %>%
    select(time_ns, hbonds, temperature) # Keep needed columns

  # Add to list
  all_hbonds_list[[temp]] <- current_hbond_data # Use temp as index
  message("  Read ", nrow(current_hbond_data), " data points for ", temp)
}

# Combine data from all temperatures
combined_hbonds <- bind_rows(all_hbonds_list)

# --- 3. Calculate Average H-bonds (Optional: for reporting) ---
avg_hbonds <- combined_hbonds %>%
  group_by(temperature) %>%
  summarise(
    avg_hbond_count = mean(hbonds, na.rm = TRUE),
    sd_hbond_count = sd(hbonds, na.rm = TRUE),
    .groups = 'drop'
  )
print("--- Average Hydrogen Bond Counts (10ns-50ns) ---")
print(avg_hbonds)

# --- 4. Generate Box Plot with Significance ---
message("\n--- Generating box plot with significance annotations ---")

# Define pairs for comparison (adjust based on expected significance)
comparison_list <- list( c("280K", "300K"), c("280K", "320K"), c("300K", "320K") )

# Calculate Y position for significance labels
y_position_for_labels <- max(combined_hbonds$hbonds, na.rm = TRUE) * 1.05 # 5% above max

# Create the box plot
hbond_boxplot_signif <- ggplot(combined_hbonds, aes(x = temperature, y = hbonds, fill = temperature)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 1.5, width = 0.6) +
  scale_fill_manual(values = temp_color_map, guide = "none") + # Use temp colors, hide legend
  # Add significance comparisons (using t-test here; consider Wilcoxon if assumptions not met)
  stat_compare_means(comparisons = comparison_list,
                     method = "t.test", # Or "wilcox.test"
                     label = "p.signif", # Show stars: *, **, ***
                     label.y = y_position_for_labels,
                     step.increase = 0.08,
                     tip.length = 0.01) +
  # Optionally add overall p-value from ANOVA/Kruskal-Wallis
  stat_compare_means(method = "anova", label.y = min(combined_hbonds$hbonds, na.rm=TRUE)*0.9) + # Adjust label.y if needed
  labs(
    title = "Distribution of Protein Hydrogen Bonds (10-50 ns)",
    # subtitle = "Significance levels: ns p > 0.05, * p <= 0.05, ** p <= 0.01, *** p <= 0.001", # Simplified levels
    x = "Temperature",
    y = "Number of H-Bonds"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) + # Expand Y axis for labels
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# --- 5. Save Plot ---
ggsave(
  filename = output_plot_filename,
  plot = hbond_boxplot_signif,
  width = 7, height = 6, dpi = 300
)
message("H-bond distribution plot saved to: ", output_plot_filename)

# --- 6. Optional: Statistical Tests (Commented out for brevity, uncomment if needed) ---
# message("\n--- Performing statistical tests (Optional) ---")
# # Check assumptions (e.g., normality, homogeneity of variances) before choosing tests
# # Example using ANOVA (if assumptions met)
# anova_result <- aov(hbonds ~ temperature, data = combined_hbonds)
# print(summary(anova_result))
# if (summary(anova_result)[[1]]$`Pr(>F)`[1] < 0.05) {
#   print(TukeyHSD(anova_result))
# }
# # Example using Kruskal-Wallis (non-parametric)
# kruskal_result <- kruskal.test(hbonds ~ temperature, data = combined_hbonds)
# print(kruskal_result)
# # if (kruskal_result$p.value < 0.05) {
# #   # Perform Dunn's test if needed (requires dunn.test package)
# #   if(require(dunn.test)){
# #        dunn_result <- dunn.test::dunn.test(combined_hbonds$hbonds, g = combined_hbonds$temperature, method="bh")
# #        print(dunn_result)
# #    } else { message("Install dunn.test package for post-hoc analysis.") }
# # }

# --- 7. End ---
message("\n--- R Script Finished ---")