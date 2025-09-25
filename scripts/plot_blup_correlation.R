# Load necessary libraries
library(dplyr)
library(ggpubr)
library(ggplot2)

# Assume 'lab_measurements' is already loaded in the environment from a previous step
# For demonstration purposes, let's create a placeholder if it doesn't exist.
if (!exists("lab_measurements")) {
  message("Creating a placeholder 'lab_measurements' data frame.")
  lab_measurements <- data.frame(
    FINNGENID = rep(paste0("FG", 1:100), each = 5),
    MEASUREMENT_VALUE_HARMONIZED = rnorm(500, mean = 100, sd = 15),
    stringsAsFactors = FALSE
  )
}

# Add columns for total measurements and median value for each FINNGENID
# and filter for individuals with more than two measurements.
lab_summary <- lab_measurements %>%
  group_by(FINNGENID) %>%
  summarise(
    n_measurements = n(),
    median_value = median(MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE)
  ) %>%
  filter(n_measurements > 2)

# Read the BLUP slopes data from the specified file
blup_slopes_file <- "00.Longitudinal_lab_values_GWAS_inputs/v.01/3001308_DF13.tsv"
if (file.exists(blup_slopes_file)) {
  blup_slopes <- read.table(blup_slopes_file, header = TRUE, sep = "\t")

  # Merge the summarized lab data with the BLUP slopes
  # The BLUP file uses 'FID' for the individual identifier, which corresponds to 'FINNGENID'
  merged_data <- inner_join(lab_summary, blup_slopes, by = c("FINNGENID" = "FID"))

  # Rename the slope column to be more descriptive.
  # The column name is based on the OMOP ID, e.g., 'X3001308_slope'
  slope_column_name <- grep("_slope$", names(merged_data), value = TRUE)[1]
  if (!is.na(slope_column_name)) {
    merged_data <- merged_data %>%
      rename(blup_slope = all_of(slope_column_name))

    # Create the scatter plot using ggpubr
    correlation_plot <- ggscatter(merged_data,
              x = "blup_slope", y = "median_value",
              add = "reg.line", conf.int = TRUE,
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "BLUP Slope",
              ylab = "Median Lab Value",
              title = "Correlation between BLUP Slopes and Median Lab Values") +
      stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
               label.x.npc = "left", label.y.npc = "top", size = 4)

    # Print the plot to the console
    print(correlation_plot)

    # Save the plot to a file
    ggsave("blup_median_correlation.png", plot = correlation_plot, width = 8, height = 6)
    message("Plot saved as blup_median_correlation.png")
  } else {
    warning("No slope column found in the BLUP data file.")
  }
} else {
  warning(paste("BLUP slopes file not found:", blup_slopes_file))
}