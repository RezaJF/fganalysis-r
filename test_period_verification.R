# Comprehensive test to verify period assignment logic
library(dplyr)
library(ggplot2)
library(ggpubr)

cat("=== TESTING PERIOD ASSIGNMENT LOGIC ===\n\n")

# Test the raw logic
cat("1. Testing time_to_drug calculation:\n")
cat("   time_to_drug = first_drug_age - EVENT_AGE\n")
cat("   Example: Drug at age 50\n")
cat("   - Event at age 45: time_to_drug = 50 - 45 = +5 (BEFORE drug)\n")
cat("   - Event at age 55: time_to_drug = 50 - 55 = -5 (AFTER drug)\n\n")

# Create test data to verify the actual behavior
test_data <- data.frame(
  FINNGENID = rep("TEST001", 6),
  EVENT_AGE = c(45, 48, 49, 51, 52, 55),
  first_drug_age = 50,
  MEASUREMENT_VALUE_HARMONIZED = c(7.0, 6.5, 6.0, 5.0, 4.5, 4.0)  # Decreasing values
)

test_data <- test_data %>%
  mutate(time_to_drug = first_drug_age - EVENT_AGE)

cat("2. Test data with time_to_drug:\n")
print(test_data)

cat("\n3. Testing period assignment with typical periods:\n")
cat("   before_period = c(-2, -0.1) (2 years to 0.1 years before drug)\n")
cat("   after_period = c(0.1, 2) (0.1 to 2 years after drug)\n\n")

# Test with what SHOULD be the correct periods based on the comment
before_period_comment <- c(-2, -0.1)  # Based on comment: negative = before
after_period_comment <- c(0.1, 2)      # Based on comment: positive = after

test_data_comment <- test_data %>%
  mutate(period_comment = case_when(
    between(time_to_drug, before_period_comment[1], before_period_comment[2]) ~ 'Before',
    between(time_to_drug, after_period_comment[1], after_period_comment[2]) ~ 'After',
    TRUE ~ NA_character_
  ))

cat("Using periods based on COMMENT (negative=before, positive=after):\n")
print(test_data_comment %>% select(EVENT_AGE, time_to_drug, MEASUREMENT_VALUE_HARMONIZED, period_comment))

# Test with what the ACTUAL logic requires
before_period_actual <- c(0.1, 5)    # Positive values for before
after_period_actual <- c(-5, -0.1)   # Negative values for after

test_data_actual <- test_data %>%
  mutate(period_actual = case_when(
    between(time_to_drug, before_period_actual[1], before_period_actual[2]) ~ 'Before',
    between(time_to_drug, after_period_actual[1], after_period_actual[2]) ~ 'After',
    TRUE ~ NA_character_
  ))

cat("\nUsing periods based on ACTUAL time_to_drug signs (positive=before, negative=after):\n")
print(test_data_actual %>% select(EVENT_AGE, time_to_drug, MEASUREMENT_VALUE_HARMONIZED, period_actual))

# Check which assignment makes sense
comment_summary <- test_data_comment %>%
  group_by(period_comment) %>%
  summarise(mean_value = mean(MEASUREMENT_VALUE_HARMONIZED), .groups = "drop") %>%
  filter(!is.na(period_comment))

actual_summary <- test_data_actual %>%
  group_by(period_actual) %>%
  summarise(mean_value = mean(MEASUREMENT_VALUE_HARMONIZED), .groups = "drop") %>%
  filter(!is.na(period_actual))

cat("\n4. Summary with COMMENT-based periods:\n")
print(comment_summary)

cat("\n5. Summary with ACTUAL-based periods:\n")
print(actual_summary)

cat("\n6. CONCLUSION:\n")
if (nrow(actual_summary) == 2) {
  before_val <- actual_summary$mean_value[actual_summary$period_actual == "Before"]
  after_val <- actual_summary$mean_value[actual_summary$period_actual == "After"]

  if (before_val > after_val) {
    cat("✓ ACTUAL periods are CORRECT: Before (", round(before_val, 2), ") > After (", round(after_val, 2), ")\n")
    cat("  This matches expected HbA1c reduction after drug treatment.\n")
  } else {
    cat("✗ ACTUAL periods seem WRONG: Before (", round(before_val, 2), ") <= After (", round(after_val, 2), ")\n")
  }
}

if (nrow(comment_summary) == 2) {
  before_val <- comment_summary$mean_value[comment_summary$period_comment == "Before"]
  after_val <- comment_summary$mean_value[comment_summary$period_comment == "After"]

  if (before_val > after_val) {
    cat("✓ COMMENT periods would be CORRECT: Before (", round(before_val, 2), ") > After (", round(after_val, 2), ")\n")
  } else {
    cat("✗ COMMENT periods are WRONG: Before (", round(before_val, 2), ") <= After (", round(after_val, 2), ")\n")
    cat("  The comment in the code is misleading!\n")
  }
}

cat("\n=== TEST COMPLETE ===\n")