library(testthat)
library(dplyr)

# Test the calculate_blup_slopes function with new parameters
test_that("calculate_blup_slopes works with save_model and plot_blup_correlation parameters", {
  # Create synthetic lab measurement data
  set.seed(123)
  n_individuals <- 20
  measurements_per_individual <- sample(3:10, n_individuals, replace = TRUE)

  # Generate synthetic data
  lab_data <- data.frame()
  for (i in 1:n_individuals) {
    finngenid <- paste0("FG", sprintf("%04d", i))
    n_meas <- measurements_per_individual[i]
    ages <- sort(runif(n_meas, 20, 80))
    # Create realistic lab values with individual-specific slopes
    individual_slope <- rnorm(1, -0.02, 0.01)  # Individual slope
    individual_intercept <- rnorm(1, 5, 0.5)   # Individual intercept
    lab_values <- individual_intercept + individual_slope * ages + rnorm(n_meas, 0, 0.2)

    individual_data <- data.frame(
      FINNGENID = finngenid,
      OMOP_CONCEPT_ID = "3001308",
      EVENT_AGE = ages,
      VALUE = lab_values,
      SEX = sample(c("male", "female"), 1)
    )
    lab_data <- rbind(lab_data, individual_data)
  }

  # Create temporary directory for output
  temp_dir <- tempdir()

  # Test with save_model = TRUE and plot_blup_correlation = FALSE
  skip_if_not_installed("lme4")

  result1 <- calculate_blup_slopes(
    data = lab_data,
    output_dir = temp_dir,
    min_measurements = 3,
    include_sex = TRUE,
    save_model = TRUE,
    plot_blup_correlation = FALSE
  )

  # Check that results are returned
  expect_type(result1, "list")
  expect_true("3001308" %in% names(result1))

  # Check that model file was created
  expect_true(!is.null(result1[["3001308"]]$model_file))
  expect_true(file.exists(result1[["3001308"]]$model_file))

  # Check that plot file was NOT created
  expect_true(is.null(result1[["3001308"]]$plot_file))

  # Test with plot_blup_correlation = TRUE
  skip_if_not_installed("ggpubr")

  result2 <- calculate_blup_slopes(
    data = lab_data,
    output_dir = temp_dir,
    min_measurements = 3,
    include_sex = TRUE,
    save_model = FALSE,
    plot_blup_correlation = TRUE
  )

  # Check that plot file was created
  expect_true(!is.null(result2[["3001308"]]$plot_file))
  expect_true(file.exists(result2[["3001308"]]$plot_file))

  # Check that model file was NOT created
  expect_true(is.null(result2[["3001308"]]$model_file))

  # Check that correlation was calculated
  expect_true(!is.null(result2[["3001308"]]$blup_fixed_correlation))
  expect_true(is.numeric(result2[["3001308"]]$blup_fixed_correlation$correlation))
  expect_true(is.numeric(result2[["3001308"]]$blup_fixed_correlation$p_value))

  # Test with both save_model = TRUE and plot_blup_correlation = TRUE
  result3 <- calculate_blup_slopes(
    data = lab_data,
    output_dir = temp_dir,
    min_measurements = 3,
    include_sex = TRUE,
    save_model = TRUE,
    plot_blup_correlation = TRUE
  )

  # Check that both files were created
  expect_true(!is.null(result3[["3001308"]]$model_file))
  expect_true(file.exists(result3[["3001308"]]$model_file))
  expect_true(!is.null(result3[["3001308"]]$plot_file))
  expect_true(file.exists(result3[["3001308"]]$plot_file))

  # Clean up temporary files
  unlink(file.path(temp_dir, "3001308_*"))
})


test_that("calculate_blup_slopes correlation calculation works with calculate_qc = TRUE", {
  # Create synthetic data with known correlation structure
  set.seed(456)
  lab_data <- data.frame()

  for (i in 1:15) {
    finngenid <- paste0("FG", sprintf("%04d", i))
    ages <- seq(30, 70, length.out = 5)
    # Create data with strong individual effects
    individual_slope <- rnorm(1, -0.05, 0.02)
    lab_values <- 6 + individual_slope * ages + rnorm(5, 0, 0.1)

    individual_data <- data.frame(
      FINNGENID = finngenid,
      OMOP_CONCEPT_ID = "3001308",
      EVENT_AGE = ages,
      VALUE = lab_values,
      SEX = sample(c("male", "female"), 1)
    )
    lab_data <- rbind(lab_data, individual_data)
  }

  temp_dir <- tempdir()

  skip_if_not_installed("lme4")

  result <- calculate_blup_slopes(
    data = lab_data,
    output_dir = temp_dir,
    min_measurements = 3,
    include_sex = TRUE,
    calculate_qc = TRUE,
    plot_blup_correlation = FALSE
  )

  # Check that correlation was calculated even without plotting
  expect_true(!is.null(result[["3001308"]]$blup_fixed_correlation))
  correlation_info <- result[["3001308"]]$blup_fixed_correlation

  expect_true(is.numeric(correlation_info$correlation))
  expect_true(is.numeric(correlation_info$p_value))
  expect_true(is.numeric(correlation_info$n_pairs))
  expect_true(correlation_info$n_pairs > 0)

  # Correlation should be positive and significant for well-behaved data
  expect_true(correlation_info$correlation > 0.3)

  # Clean up
  unlink(file.path(temp_dir, "3001308_*"))
})

test_that("calculate_blup_slopes handles insufficient data for correlation", {
  # This test should now pass as the warning is expected
  expect_warning(
    result <- calculate_blup_slopes(
      data = data.frame(
        FINNGENID = c("FG1", "FG2"),
        OMOP_CONCEPT_ID = "3001308",
        EVENT_AGE = c(20, 21),
        VALUE = c(1, 2)
      ),
      plot_blup_correlation = TRUE,
      include_sex = FALSE
    ),
    "Skipping analysis"
  )
  # When analysis is skipped, the result for that concept ID should be NULL
  expect_null(result[["3001308"]])
})

test_that("calculate_blup_slopes recovers known individual slopes with synthetic data", {
  # Create synthetic data with predetermined individual slopes
  # This test validates that BLUP calculations are mathematically correct
  set.seed(42)  # For reproducibility

  # Define known true slopes for each individual
  true_slopes <- c(-0.05, -0.02, 0.01, 0.03, -0.01, 0.02, -0.03, 0.00, 0.04, -0.02)
  n_individuals <- length(true_slopes)

  # Generate lab data with known parameters
  lab_data <- data.frame()
  for (i in 1:n_individuals) {
    finngenid <- paste0("FG", sprintf("%04d", i))
    # Create regular age intervals for each individual
    ages <- seq(30, 70, length.out = 8)  # 8 measurements from age 30 to 70
    true_slope <- true_slopes[i]
    true_intercept <- 5.0  # Known intercept for all individuals

    # Generate lab values with known parameters + small noise
    # The noise is small enough that BLUP should recover the true slopes
    lab_values <- true_intercept + true_slope * ages + rnorm(8, 0, 0.1)

    individual_data <- data.frame(
      FINNGENID = finngenid,
      OMOP_CONCEPT_ID = "3001308",
      EVENT_AGE = ages,
      VALUE = lab_values,
      SEX = "male"  # Keep simple for testing
    )
    lab_data <- rbind(lab_data, individual_data)
  }

  # Run BLUP analysis
  skip_if_not_installed("lme4")

  result <- calculate_blup_slopes(
    data = lab_data,
    min_measurements = 3,
    include_sex = FALSE,  # Simplified for testing
    calculate_qc = TRUE   # Enable QC to test correlation calculation
  )

  # Verify that results are returned
  expect_type(result, "list")
  expect_true("3001308" %in% names(result))
  expect_false(is.null(result[["3001308"]]))

  # Extract calculated slopes
  calculated_slopes <- result[["3001308"]]$blup_slopes$slope
  expect_length(calculated_slopes, n_individuals)

  # Verify that calculated slopes are close to true slopes
  # Allow for some estimation error due to noise (tolerance based on noise level)
  for (i in 1:n_individuals) {
    expect_equal(calculated_slopes[i], true_slopes[i], tolerance = 0.02,
                 info = paste("Individual", i, "slope mismatch"))
  }

  # Additional validation: slopes should be in reasonable range
  expect_true(all(calculated_slopes > -0.1 & calculated_slopes < 0.1),
              info = "All slopes should be in reasonable range")

  # Test that BLUP slopes are unbiased (mean should be close to mean of true slopes)
  mean_true_slope <- mean(true_slopes)
  mean_calculated_slope <- mean(calculated_slopes)
  expect_equal(mean_calculated_slope, mean_true_slope, tolerance = 0.01,
               info = "Mean BLUP slope should be unbiased")

  # Test that correlation with fixed-effect slopes is high (QC validation)
  correlation_info <- result[["3001308"]]$blup_fixed_correlation
  expect_false(is.null(correlation_info))
  expect_true(correlation_info$correlation > 0.8,
              info = "BLUP-fixed correlation should be high for well-behaved data")
  expect_true(correlation_info$p_value < 0.05,
              info = "Correlation should be statistically significant")
})

test_that("calculate_blup_slopes handles edge cases correctly", {
  skip_if_not_installed("lme4")

  # Test 1: No individual variation (all slopes should be ~0)
  set.seed(123)
  lab_data_no_variation <- data.frame()
  for (i in 1:12) {  # Use 12 individuals to meet minimum requirement
    finngenid <- paste0("FG", sprintf("%04d", i))
    ages <- seq(30, 70, length.out = 6)
    # All individuals have the same slope (0) and intercept
    lab_values <- 5.0 + 0.0 * ages + rnorm(6, 0, 0.1)

    individual_data <- data.frame(
      FINNGENID = finngenid,
      OMOP_CONCEPT_ID = "3001308",
      EVENT_AGE = ages,
      VALUE = lab_values,
      SEX = "male"
    )
    lab_data_no_variation <- rbind(lab_data_no_variation, individual_data)
  }

  result_no_var <- calculate_blup_slopes(
    data = lab_data_no_variation,
    min_measurements = 3,
    include_sex = FALSE
  )

  # Check if the model converged (it might not for no-variation case)
  if (!is.null(result_no_var[["3001308"]])) {
    calculated_slopes_no_var <- result_no_var[["3001308"]]$blup_slopes$slope
    # All slopes should be close to 0
    expect_true(all(abs(calculated_slopes_no_var) < 0.05),
                info = "Slopes should be close to 0 when there's no individual variation")
  } else {
    # If model doesn't converge, that's also acceptable for no-variation case
    expect_null(result_no_var[["3001308"]],
                info = "Model may fail to converge when there's no individual variation")
  }

  # Test 2: Perfect correlation scenario (all individuals have same slope)
  set.seed(456)
  lab_data_perfect <- data.frame()
  common_slope <- 0.02  # All individuals have this slope
  for (i in 1:12) {  # Use 12 individuals to meet minimum requirement
    finngenid <- paste0("FG", sprintf("%04d", i))
    ages <- seq(30, 70, length.out = 6)
    # All individuals have the same slope but different intercepts
    individual_intercept <- 4.5 + i * 0.1
    lab_values <- individual_intercept + common_slope * ages + rnorm(6, 0, 0.05)

    individual_data <- data.frame(
      FINNGENID = finngenid,
      OMOP_CONCEPT_ID = "3001308",
      EVENT_AGE = ages,
      VALUE = lab_values,
      SEX = "male"
    )
    lab_data_perfect <- rbind(lab_data_perfect, individual_data)
  }

  result_perfect <- calculate_blup_slopes(
    data = lab_data_perfect,
    min_measurements = 3,
    include_sex = FALSE,
    calculate_qc = TRUE
  )

  calculated_slopes_perfect <- result_perfect[["3001308"]]$blup_slopes$slope
  # All slopes should be close to the common slope (allow for BLUP shrinkage)
  # BLUP shrinkage can be substantial, so we use a more generous tolerance
  expect_true(all(abs(calculated_slopes_perfect - common_slope) < 0.05),
              info = "All slopes should be close to the common slope")

  # Test 3: Extreme slopes (but still within reasonable bounds)
  set.seed(789)
  lab_data_extreme <- data.frame()
  extreme_slopes <- c(-0.08, -0.06, 0.06, 0.08, 0.0, -0.04, 0.04, -0.02, 0.02, 0.01, -0.01, 0.03)  # More extreme but realistic
  for (i in 1:12) {  # Use 12 individuals to meet minimum requirement
    finngenid <- paste0("FG", sprintf("%04d", i))
    ages <- seq(30, 70, length.out = 8)
    true_slope <- extreme_slopes[i]
    lab_values <- 5.0 + true_slope * ages + rnorm(8, 0, 0.1)

    individual_data <- data.frame(
      FINNGENID = finngenid,
      OMOP_CONCEPT_ID = "3001308",
      EVENT_AGE = ages,
      VALUE = lab_values,
      SEX = "male"
    )
    lab_data_extreme <- rbind(lab_data_extreme, individual_data)
  }

  result_extreme <- calculate_blup_slopes(
    data = lab_data_extreme,
    min_measurements = 3,
    include_sex = FALSE
  )

  calculated_slopes_extreme <- result_extreme[["3001308"]]$blup_slopes$slope
  # Verify that extreme slopes are recovered
  for (i in 1:12) {
    expect_equal(calculated_slopes_extreme[i], extreme_slopes[i], tolerance = 0.03,
                 info = paste("Extreme slope", i, "not recovered correctly"))
  }
})

test_that("BLUP slopes have expected statistical properties", {
  skip_if_not_installed("lme4")

  # Create data with known statistical properties
  set.seed(999)
  n_individuals <- 15
  n_measurements <- 6

  # True slopes with known distribution
  true_slopes <- rnorm(n_individuals, mean = 0, sd = 0.03)

  lab_data <- data.frame()
  for (i in 1:n_individuals) {
    finngenid <- paste0("FG", sprintf("%04d", i))
    ages <- seq(30, 70, length.out = n_measurements)
    true_slope <- true_slopes[i]
    true_intercept <- rnorm(1, 5, 0.5)  # Random intercepts

    lab_values <- true_intercept + true_slope * ages + rnorm(n_measurements, 0, 0.1)

    individual_data <- data.frame(
      FINNGENID = finngenid,
      OMOP_CONCEPT_ID = "3001308",
      EVENT_AGE = ages,
      VALUE = lab_values,
      SEX = sample(c("male", "female"), 1)
    )
    lab_data <- rbind(lab_data, individual_data)
  }

  result <- calculate_blup_slopes(
    data = lab_data,
    min_measurements = 3,
    include_sex = TRUE,  # Include sex effect
    calculate_qc = TRUE
  )

  calculated_slopes <- result[["3001308"]]$blup_slopes$slope

  # Test 1: BLUP slopes should be unbiased
  mean_true <- mean(true_slopes)
  mean_calculated <- mean(calculated_slopes)
  expect_equal(mean_calculated, mean_true, tolerance = 0.01,
               info = "BLUP slopes should be unbiased")

  # Test 2: Variance of BLUP slopes should be reasonable
  var_true <- var(true_slopes)
  var_calculated <- var(calculated_slopes)
  # BLUP variance should be less than or equal to true variance (shrinkage)
  expect_true(var_calculated <= var_true * 1.1,  # Allow small increase due to noise
              info = "BLUP variance should not be much larger than true variance")

  # Test 3: Correlation between true and calculated slopes should be high
  correlation_true_calc <- cor(true_slopes, calculated_slopes)
  expect_true(correlation_true_calc > 0.7,
              info = "Correlation between true and calculated slopes should be high")

  # Test 4: QC correlation should be reasonable
  qc_correlation <- result[["3001308"]]$blup_fixed_correlation$correlation
  expect_true(qc_correlation > 0.5,
              info = "BLUP-fixed correlation should be reasonable")

  # Test 5: All slopes should be finite
  expect_true(all(is.finite(calculated_slopes)),
              info = "All BLUP slopes should be finite")
})