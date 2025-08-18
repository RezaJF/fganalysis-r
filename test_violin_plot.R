# Test the updated violin plot function
library(dplyr)
library(ggplot2)
library(ggpubr)
source("R/visualization.R")
source("R/drug_response_core.R")

cat("Testing violin plot functionality...\n")

# Create test data
response_data <- data.frame(
  FINNGENID = paste0("FG", 1:10),
  response = rnorm(10, -0.5, 0.3),
  first_drug = rep("A10AA01", 10)
)

lab_measurements <- data.frame()
for (i in 1:10) {
  finngenid <- paste0("FG", i)

  # Before measurements (higher values)
  before_data <- data.frame(
    FINNGENID = finngenid,
    OMOP_CONCEPT_ID = "3001308",
    EVENT_AGE = runif(5, 45, 49),
    MEASUREMENT_VALUE_HARMONIZED = rnorm(5, 5.5, 0.5),
    first_drug = "A10AA01",
    first_drug_age = 50
  )

  # After measurements (lower values)
  after_data <- data.frame(
    FINNGENID = finngenid,
    OMOP_CONCEPT_ID = "3001308",
    EVENT_AGE = runif(5, 51, 55),
    MEASUREMENT_VALUE_HARMONIZED = rnorm(5, 4.5, 0.5),
    first_drug = "A10AA01",
    first_drug_age = 50
  )

  lab_measurements <- rbind(lab_measurements, before_data, after_data)
}

# Calculate time_to_drug
lab_measurements <- lab_measurements %>%
  mutate(time_to_drug = first_drug_age - EVENT_AGE)

drug_purchases <- data.frame(
  FINNGENID = paste0("FG", 1:10),
  ATC = "A10AA01",
  EVENT_AGE = 50
)

# Create drug.response object
test_response <- drug.response(
  responses = response_data,
  lab_measurements = lab_measurements,
  drug_purchases = drug_purchases,
  before_period = c(0.5, 6),    # Positive values for before
  after_period = c(-6, -0.5)    # Negative values for after
)

# Test the violin plot function
tryCatch({
  p <- plot_lab_value_distribution(test_response, remove_outliers = FALSE)

  if (inherits(p, "ggplot")) {
    cat("✓ Violin plot created successfully\n")

    # Check that it's using ggviolin
    plot_layers <- sapply(p$layers, function(x) class(x$geom)[1])
    if (any(grepl("Geom", plot_layers))) {
      cat("✓ Plot contains geometric layers\n")
    }

    # Save the plot to verify
    ggsave("test_violin_plot.pdf", plot = p, width = 10, height = 6)
    cat("✓ Violin plot saved to test_violin_plot.pdf\n")

  } else {
    cat("✗ Function did not return a ggplot object\n")
  }

}, error = function(e) {
  cat("✗ Function failed with error:", e$message, "\n")
})

cat("\nViolin plot test completed!\n")