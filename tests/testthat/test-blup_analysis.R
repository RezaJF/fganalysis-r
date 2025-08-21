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
      MEASUREMENT_VALUE_HARMONIZED = lab_values,
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
      MEASUREMENT_VALUE_HARMONIZED = lab_values,
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
  # Create data with very few individuals
  lab_data <- data.frame(
    FINNGENID = rep(c("FG001", "FG002"), each = 3),
    OMOP_CONCEPT_ID = "3001308",
    EVENT_AGE = c(25, 35, 45, 30, 40, 50),
    MEASUREMENT_VALUE_HARMONIZED = rnorm(6, 5, 0.5),
    SEX = rep("male", 6)
  )

  temp_dir <- tempdir()

  skip_if_not_installed("lme4")

  # Should warn about insufficient data but still run
  expect_warning(
    result <- calculate_blup_slopes(
      data = lab_data,
      output_dir = temp_dir,
      min_measurements = 3,
      include_sex = TRUE,
      calculate_qc = TRUE
    ),
    "Insufficient data for correlation analysis"
  )

  # Should still return results but without correlation info
  expect_true("3001308" %in% names(result))
  expect_true(is.null(result[["3001308"]]$blup_fixed_correlation))

  # Clean up
  unlink(file.path(temp_dir, "3001308_*"))
})