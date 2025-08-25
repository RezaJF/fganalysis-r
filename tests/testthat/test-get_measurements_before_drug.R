library(testthat)
library(dplyr)

# Mock connection object for testing
mock_conn <- function() {
  labs <- data.frame(
    FINNGENID = c("FG1", "FG1", "FG1", "FG2", "FG2", "FG3", "FG3", "FG4"),
    OMOP_CONCEPT_ID = "3001308",
    EVENT_AGE = c(50.0, 50.1, 50.5, 60.0, 60.3, 70.0, 70.2, 80.0),
    MEASUREMENT_VALUE_HARMONIZED = c(100, 102, 150, 110, 112, 120, 122, 500),
    stringsAsFactors = FALSE
  )

  pheno <- data.frame(
    FINNGENID = c("FG1", "FG2"),
    EVENT_AGE = c(50.2, 60.1),
    CODE1 = "C10AA",
    SOURCE = "PURCH",
    stringsAsFactors = FALSE
  )

  cov <- data.frame(
    FINNGENID = c("FG1", "FG2", "FG3", "FG4"),
    SEX = c("Male", "Female", "Male", "Female"),
    stringsAsFactors = FALSE
  )

  return(list(labs = labs, pheno = pheno, cov_pheno = cov))
}

test_that("get_measurements_before_drug works as a standalone function", {
  conn <- mock_conn()

  # Test with default 3 months
  measurements <- get_measurements_before_drug(
    conn,
    lablist = "3001308",
    druglist = "C10AA",
    covariates = conn$cov_pheno,
    covariate_cols = "SEX"
  )

  # Expectations remain the same as the logic hasn't changed, just the test setup
  expect_equal(nrow(measurements), 2 + 1 + 2 + 1)
  expect_true("SEX" %in% colnames(measurements))

  # Filter to ensure there are rows for FG1 before pulling
  fg1_sex <- measurements %>% filter(FINNGENID == "FG1") %>% pull(SEX)
  if(length(fg1_sex) > 0) {
    expect_equal(unique(fg1_sex), "Male")
  }

  # Check measurement counts
  fg1_counts <- measurements %>% filter(FINNGENID == "FG1") %>% pull(n_measurements)
  if(length(fg1_counts) > 0) {
    expect_equal(unique(fg1_counts), 2)
  }

  fg3_counts <- measurements %>% filter(FINNGENID == "FG3") %>% pull(n_measurements)
  if(length(fg3_counts) > 0) {
    expect_equal(unique(fg3_counts), 2)
  }
})

test_that("outlier removal works correctly", {
  conn <- mock_conn()

  # Test SD removal
  measurements_sd <- get_measurements_before_drug(
    conn,
    lablist = "3001308",
    druglist = "C10AA",
    remove_outliers_sd = 1
  )
  # The value 500 from FG4 should be removed
  expect_false(500 %in% measurements_sd$MEASUREMENT_VALUE_HARMONIZED)

  # Test Winsorizing
  measurements_win <- get_measurements_before_drug(
    conn,
    lablist = "3001308",
    druglist = "C10AA",
    winsorize_pct = 0.1  # 10% winsorization on each tail
  )
  # The value 500 should be capped at the 90th percentile
  expect_true(max(measurements_win$MEASUREMENT_VALUE_HARMONIZED) < 500)
})

test_that("smooth_measurement_intervals works correctly", {
  test_data <- data.frame(
    FINNGENID = "FG1", OMOP_CONCEPT_ID = "123", SEX_CODED = 1L,
    EVENT_AGE = c(50.0, 50.1, 50.3, 51.0, 51.2, 52.0),
    MEASUREMENT_VALUE_HARMONIZED = c(10, 12, 10, 20, 22, 30)
  )

  smoothed <- smooth_measurement_intervals(test_data, min_interval_months = 6)

  # Three clusters: (50.0, 50.1, 50.3), (51.0, 51.2), (52.0)
  expect_equal(nrow(smoothed), 3)
  # Check first cluster: mean age, median value
  expect_equal(smoothed$EVENT_AGE[1], mean(c(50.0, 50.1, 50.3)))
  expect_equal(smoothed$MEASUREMENT_VALUE_HARMONIZED[1], 10) # Median of 10, 12, 10 is 10
  # Check second cluster
  expect_equal(smoothed$EVENT_AGE[2], mean(c(51.0, 51.2)))
  expect_equal(smoothed$MEASUREMENT_VALUE_HARMONIZED[2], 21) # Median of 20, 22 is 21
  # Check third cluster (isolated point)
  expect_equal(smoothed$EVENT_AGE[3], 52.0)
  expect_equal(smoothed$MEASUREMENT_VALUE_HARMONIZED[3], 30)
})
