# Test the fixed period logic
library(dplyr)
library(ggplot2)
library(ggpubr)
source("R/visualization.R")
source("R/drug_response_core.R")

cat("Testing FIXED period logic...\n\n")

# Create test data with clear before/after pattern
response_data <- data.frame(
  FINNGENID = paste0("FG", 1:5),
  response = rep(-1.5, 5),  # Negative response (after < before)
  first_drug = rep("A10AA01", 5)
)

# Create lab measurements
lab_measurements <- data.frame()
for (i in 1:5) {
  finngenid <- paste0("FG", i)

  # Before measurements (6 months to 1 year before drug, higher HbA1c values)
  before_data <- data.frame(
    FINNGENID = finngenid,
    OMOP_CONCEPT_ID = "3001308",
    EVENT_AGE = runif(3, 49, 49.5),  # 6-12 months before age 50
    MEASUREMENT_VALUE_HARMONIZED = rnorm(3, 7.5, 0.3),  # High HbA1c
    first_drug = "A10AA01",
    first_drug_age = 50
  )

  # After measurements (1 month to 1 year after drug, lower HbA1c values)
  after_data <- data.frame(
    FINNGENID = finngenid,
    OMOP_CONCEPT_ID = "3001308",
    EVENT_AGE = runif(3, 50.1, 51),  # 1 month to 1 year after age 50
    MEASUREMENT_VALUE_HARMONIZED = rnorm(3, 6.0, 0.3),  # Lower HbA1c
    first_drug = "A10AA01",
    first_drug_age = 50
  )

  lab_measurements <- rbind(lab_measurements, before_data, after_data)
}

# Calculate time_to_drug
lab_measurements <- lab_measurements %>%
  mutate(time_to_drug = first_drug_age - EVENT_AGE)

drug_purchases <- data.frame(
  FINNGENID = paste0("FG", 1:5),
  ATC = "A10AA01",
  EVENT_AGE = 50
)

# Create drug.response object with typical periods
# Using the convention that the function expects
test_response <- drug.response(
  responses = response_data,
  lab_measurements = lab_measurements,
  drug_purchases = drug_purchases,
  before_period = c(-1, 0),    # 1 year to 0 before drug (negative values)
  after_period = c(1/12, 1)    # 1 month to 1 year after drug (positive values)
)

cat("Period definitions used:\n")
cat("  before_period = c(-1, 0) (convention: 1 year to 0 before drug)\n")
cat("  after_period = c(1/12, 1) (convention: 1 month to 1 year after drug)\n\n")

# Test the function
tryCatch({
  p <- plot_lab_value_distribution(test_response, remove_outliers = FALSE)

  if (inherits(p, "ggplot")) {
    cat("✓ Function runs without errors\n")

    # Check the data used for plotting
    plot_data <- p$data

    cat("\nPeriod assignments in plot data:\n")
    period_summary <- plot_data %>%
      group_by(period) %>%
      summarise(
        n = n(),
        mean_value = mean(MEASUREMENT_VALUE_HARMONIZED),
        mean_time_to_drug = mean(time_to_drug),
        mean_event_age = mean(EVENT_AGE),
        .groups = "drop"
      )
    print(period_summary)

    # Verify logic: Before should have higher HbA1c values than After
    before_mean <- period_summary$mean_value[period_summary$period == "Before"]
    after_mean <- period_summary$mean_value[period_summary$period == "After"]

    cat("\n")
    if (!is.na(before_mean) && !is.na(after_mean)) {
      if (before_mean > after_mean) {
        cat("✓ CORRECT: Before HbA1c (", round(before_mean, 2), ") > After HbA1c (", round(after_mean, 2), ")\n")
        cat("  This matches expected HbA1c reduction after diabetes drug treatment!\n")
      } else {
        cat("✗ ERROR: Before HbA1c (", round(before_mean, 2), ") <= After HbA1c (", round(after_mean, 2), ")\n")
        cat("  This is opposite of what we expect - periods may still be mixed up!\n")
      }
    } else {
      cat("✗ ERROR: Could not find both Before and After periods in the data\n")
    }

    # Save the plot
    ggsave("test_fixed_periods_plot.pdf", plot = p, width = 10, height = 6)
    cat("\n✓ Plot saved to test_fixed_periods_plot.pdf\n")

  } else {
    cat("✗ Function did not return a ggplot object\n")
  }

}, error = function(e) {
  cat("✗ Function failed with error:", e$message, "\n")
})

cat("\nFixed period test completed!\n")