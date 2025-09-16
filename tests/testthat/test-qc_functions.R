library(testthat)
library(fganalysis)

test_that("quantile_normalize handles various inputs correctly", {
  # Test with normal data
  set.seed(123)
  x <- rnorm(100, mean = 10, sd = 2)
  x_norm <- quantile_normalize(x)

  expect_equal(length(x_norm), length(x))
  expect_true(abs(mean(x_norm, na.rm = TRUE)) < 0.1)  # Should be close to 0
  expect_true(abs(sd(x_norm, na.rm = TRUE) - 1) < 0.2)  # Should be close to 1

  # Test with NA values
  x_na <- c(1, 2, NA, 4, 5, NA)
  x_na_norm <- quantile_normalize(x_na)

  expect_equal(which(is.na(x_na)), which(is.na(x_na_norm)))
  expect_equal(length(x_na_norm), length(x_na))

  # Test with all NA
  x_all_na <- rep(NA, 5)
  expect_equal(quantile_normalize(x_all_na), x_all_na)

  # Test with empty vector
  expect_equal(quantile_normalize(numeric(0)), numeric(0))

  # Test with single value
  expect_equal(quantile_normalize(5), 0)
})

test_that("calculate_fixed_slopes works correctly", {
  # Create test data
  set.seed(123)
  test_data <- data.frame(
    FINNGENID = rep(c("ID1", "ID2", "ID3"), each = 10),
    EVENT_AGE = rep(1:10, 3),
    MEASUREMENT_VALUE_HARMONIZED = c(
      1:10 + rnorm(10, 0, 0.5),  # Positive slope
      20 - 1:10 + rnorm(10, 0, 0.5),  # Negative slope
      rep(15, 10) + rnorm(10, 0, 0.5)  # No slope
    )
  )

  # Calculate slopes
  slopes <- calculate_fixed_slopes(test_data, min_measurements = 5)

  expect_equal(nrow(slopes), 3)
  expect_true("FINNGENID" %in% colnames(slopes))
  expect_true("fixed_slope" %in% colnames(slopes))

  # Check slope directions
  expect_true(slopes$fixed_slope[slopes$FINNGENID == "ID1"] > 0)
  expect_true(slopes$fixed_slope[slopes$FINNGENID == "ID2"] < 0)
  expect_true(abs(slopes$fixed_slope[slopes$FINNGENID == "ID3"]) < 0.5)

  # Test with insufficient measurements
  test_data_few <- test_data[test_data$FINNGENID == "ID1",][1:3,]
  slopes_few <- calculate_fixed_slopes(test_data_few, min_measurements = 5)
  expect_equal(nrow(slopes_few), 0)
})

test_that("process_variance_files handles file operations correctly", {
  # Create temporary directory and test files
  temp_dir <- tempdir()

  # Create test variance file
  test_variance_data <- data.frame(
    FID = paste0("ID", 1:100),
    IID = paste0("ID", 1:100),
    test_variance = c(rnorm(90, mean = 5, sd = 2), rep(NA, 10))
  )

  test_file <- file.path(temp_dir, "test_variance.tsv")
  write.table(test_variance_data, file = test_file, sep = "\t",
              row.names = FALSE, quote = FALSE)

  # Process files without plots
  summary_table <- process_variance_files(temp_dir, generate_plots = FALSE,
                                          save_normalized = TRUE)

  expect_true(!is.null(summary_table))
  expect_equal(nrow(summary_table), 1)
  expect_true(all(c("File", "Total_records", "Non_missing", "Orig_SD",
                    "Norm_SD") %in% colnames(summary_table)))

  # Check that normalized file was created
  norm_file <- file.path(temp_dir, "test_variance_qnorm.tsv")
  expect_true(file.exists(norm_file))

  # Read normalized file and check
  norm_data <- read.table(norm_file, header = TRUE, sep = "\t")
  expect_true("test_variance_qnorm" %in% colnames(norm_data))

  # Clean up
  unlink(test_file)
  unlink(norm_file)
})

test_that("QC correlation calculation works in calculate_blup_slopes", {
  skip("Requires full drug_response object setup")

  # This test would require a complete mock setup of drug_response object
  # which is complex. In practice, this would be tested with real data
  # or a comprehensive mock object.
})

test_that("calculate_blup_slopes works with direct lab measurement input", {
  # Create test lab measurement data
  set.seed(123)
  n_individuals <- 20
  measurements_per_person <- 15

  # Generate test data
  test_lab_data <- expand.grid(
    FINNGENID = paste0("ID", 1:n_individuals),
    measurement = 1:measurements_per_person
  )

  # Add required columns
  test_lab_data$OMOP_CONCEPT_ID <- "3001308"  # LDL cholesterol
  test_lab_data$EVENT_AGE <- test_lab_data$measurement +
    rep(runif(n_individuals, 40, 60), each = measurements_per_person)

  # Create different slopes for different individuals
  individual_slopes <- runif(n_individuals, -0.5, 0.5)
  individual_intercepts <- runif(n_individuals, 2, 5)

  test_lab_data$MEASUREMENT_VALUE_HARMONIZED <- NA
  for (i in 1:n_individuals) {
    idx <- test_lab_data$FINNGENID == paste0("ID", i)
    test_lab_data$MEASUREMENT_VALUE_HARMONIZED[idx] <-
      individual_intercepts[i] +
      individual_slopes[i] * test_lab_data$EVENT_AGE[idx] +
      rnorm(sum(idx), 0, 0.2)
  }

  # Add SEX column
  test_lab_data$SEX <- rep(sample(c("M", "F"), n_individuals, replace = TRUE),
                           each = measurements_per_person)

  # Create temporary output directory
  temp_dir <- tempdir()

  # Test basic functionality without drug response
  expect_no_error({
    blup_results <- calculate_blup_slopes(
      data = test_lab_data,
      output_dir = temp_dir,
      min_measurements = 5,
      include_sex = TRUE,
      calculate_qc = TRUE  # Test QC with direct input
    )
  })

  # Check that results were created
  expect_true(length(blup_results) > 0)
  expect_true("3001308" %in% names(blup_results))

    # Check output file exists
  output_file <- file.path(temp_dir, "3001308_DF13_blup.tsv")
  expect_true(file.exists(output_file))

  # Read and check output
  output_data <- read.table(output_file, header = TRUE, sep = "\t")
  expect_equal(nrow(output_data), n_individuals)
  expect_true("X3001308_slope" %in% colnames(output_data))
  expect_true("X3001308_fixed_slope" %in% colnames(output_data))  # QC was enabled

  # Test without SEX column
  test_lab_data_no_sex <- test_lab_data[, !colnames(test_lab_data) %in% "SEX"]

  expect_no_error({
    blup_results_no_sex <- calculate_blup_slopes(
      data = test_lab_data_no_sex,
      output_dir = temp_dir,
      min_measurements = 5,
      include_sex = FALSE  # Must set to FALSE when no SEX column
    )
  })

  # Test that the function stops if SEX column is missing and include_sex is TRUE
  # The error message should now check for either SEX or SEX_IMPUTED
  expect_error(
    calculate_blup_slopes(test_lab_data_no_sex, include_sex = TRUE),
    "SEX or SEX_IMPUTED column not found"
  )

  # Test warning for unsupported options
  expect_warning({
    calculate_blup_slopes(
      data = test_lab_data,
      output_dir = temp_dir,
      drug_exposed_only = TRUE,  # Not supported with direct input
      calculate_post_variance = TRUE  # Not supported with direct input
    )
  }, "not supported with direct lab measurement input")

  # Clean up
  unlink(list.files(temp_dir, pattern = "_DF13\\.tsv$", full.names = TRUE))
})