library(testthat)


# Load the functions from the package
#source("R/drug_response_functions.R")

# Test the drug.response function
test_that("drug.response creates the correct object", {
  response <- data.frame(FINNGENID = c(1, 2), response = c(1, 2))
  lab_measurements <- data.frame(FINNGENID = c(1, 2), MEASUREMENT_VALUE_HARMONIZED = c(10, 20))
  drug_purchases <- data.frame(FINNGENID = c(1, 2), ATC = c("A01", "A02"))

  result <- drug.response(response, lab_measurements, drug_purchases, c(-1, -0.5), c(0.5, 1))

  expect_s3_class(result, "drug.reponse")
  expect_equal(result$response, response)
  expect_equal(result$all_measurements, lab_measurements)
  expect_equal(result$all_drug_purchases, drug_purchases)
})

# Test the generate_response_summary function
test_that("generate_response_summary calculates correct summaries", {
lab_measurements <- data.frame(FINNGENID = c("FG1", "FG1", "FG1", "FG1", "FG1", "FG2", "FG2", "FG2", "FG2", "FG3", "FG3"),
                                  EVENT_AGE = c(21.1, 20, 20.5, 21.5, 22.0, 34, 34.4, 33.5, 35.0, 40, 40.5),
                                  MEASUREMENT_VALUE_HARMONIZED = c(10, 20, 42, 15, 12, 30, 44, 25, 50, 120, 38),
                                  first_drug = c("A01", "A01", "A01", "A01", "A01", "A02", "A02", "A02", "A02", "A03", "A03"),
                                  first_drug_age = c(21.05, 21.05, 21.05, 21.05, 21.05, 34.2, 34.2, 34.2, 34.2, 35, 35))
  lab_measurements <-  lab_measurements %>% mutate(time_to_drug = first_drug_age - EVENT_AGE)

  before_period <- c(-1.5, 0)
  after_period <- c(0.00001, 1.5)

  result <- generate_response_summary(lab_measurements, before_period, after_period)

  expect_equal(nrow(result), 2)
  expect_equal(result$before, c(31, 27.5))
  expect_equal(result$after, c(12, 47))
  expect_equal(result$response, c(-19, 19.5))
})

# Test the quant_text function
test_that("quant_text formats quantiles correctly", {
  vector <- c(1, 2, 3, 4, 5)
  result <- quant_text(vector)

  expect_true(grepl("0%:", result))
  expect_true(grepl("100%:", result))
})

# Test the create_drug_response function
test_that("create_drug_response returns the correct structure", {
  kanta <- data.frame(
    FINNGENID = c("FG1", "FG1", "FG1", "FG1", "FG2", "FG2", "FG2", "FG2","FG3", "FG3"),
    OMOP_CONCEPT_ID = c("lab1", "lab1", "lab1", "lab1", "lab2", "lab2", "lab2", "lab2", "lab2", "lab2"),
    EVENT_AGE = c(20.6, 20.7, 20.8, 21.5, 19.5, 19.6, 19.7, 20.5,25, 25.5),
    MEASUREMENT_VALUE_HARMONIZED = c(15, 16, 17, 25, 8, 9, 10, 40, 50, 38))

  phenos <- data.frame(
    FINNGENID = c("FG1", "FG2","FG3"),
    SOURCE = c("PURCH", "PURCH","PURCH"),
    APPROX_EVENT_DAY = as.Date(c("2015-07-17" , "2015-07-18", "2015-07-19")),
    CODE1 = c("A01", "A02", "A02"),
    CODE2 = c("", "", ""),
    CODE3 = c("", "", ""),
    CODE4 = c("1", "1", "1"),
    EVENT_AGE = c(21.0, 20.0, 35)
  )

  conn <-   fg_data_connection(list(pheno = phenos, labs = kanta))

  lablist <- c("lab1", "lab2")
  druglist <- c("A01", "A02")

  result <- create_drug_response(conn, lablist, druglist, c(-1, 0), c(0.1, 1))

  expect_s3_class(result, "drug.reponse")
  expect_equal(nrow(result$response), 2)
  expect_equal(result$response$FINNGENID, c("FG1", "FG2"))
  expect_equal(result$response$before, c(16, 9))
  expect_equal(result$response$after, c(25, 40))
  expect_equal(result$response$response, c(9, 31))
})

test_that("create_drug_response removes outliers correctly", {
  kanta <- data.frame(
    FINNGENID = c("FG1", "FG1", "FG1", "FG1", "FG2", "FG2", "FG2", "FG2", "FG3", "FG3", "FG4", "FG4"),
    OMOP_CONCEPT_ID = c("lab1", "lab1", "lab1", "lab1", "lab1", "lab1", "lab1", "lab1", "lab1", "lab1", "lab1", "lab1"),
    EVENT_AGE = c(20.6, 20.7, 20.8, 21.5, 19.5, 19.6, 19.7, 20.5, 25, 25.5, 30, 30.5),
    MEASUREMENT_VALUE_HARMONIZED = c(15, 16, 17, 25, 8, 9, 10, 40, 50, 38, 1000, 1001) # Add a clear outlier
  )

  phenos <- data.frame(
    FINNGENID = c("FG1", "FG2", "FG3", "FG4"),
    SOURCE = c("PURCH", "PURCH", "PURCH", "PURCH"),
    CODE1 = c("A01", "A02", "A02", "A01"),
    CODE2 = c("", "", "", ""),
    CODE3 = c("", "", "", ""),
    CODE4 = c("1", "1", "1", "1"),
    EVENT_AGE = c(21.0, 20.0, 35.0, 31.0)
  )

  lablist <- c("lab1")
  druglist <- c("A01", "A02")

  # Test with outlier removal
  result <- create_drug_response(kanta, phenos, lablist, druglist, c(-1, 0), c(0.1, 1), remove_outliers_sd = 1)

  # FG4 has the outlier and should be removed from the response
  expect_s3_class(result, "drug.reponse")
  expect_equal(nrow(result$responses), 2)
  expect_false("FG4" %in% result$responses$FINNGENID)

  # Test that it throws an error with invalid input
  expect_error(create_drug_response(kanta, phenos, lablist, druglist, c(-1, 0), c(0.1, 1), remove_outliers_sd = 7))
  expect_error(create_drug_response(kanta, phenos, lablist, druglist, c(-1, 0), c(0.1, 1), remove_outliers_sd = "a"))
})

test_that("summarize_drug_purchases_upset creates a plot file", {
  # Create a dummy drug.response object
  response <- data.frame(FINNGENID = c(1, 2), response = c(1, 2))
  lab_measurements <- data.frame(FINNGENID = c(1, 2), MEASUREMENT_VALUE_HARMONIZED = c(10, 20))
  drug_purchases <- data.frame(
    FINNGENID = c("FG1", "FG1", "FG2", "FG2", "FG3"),
    ATC = c("A01", "A02", "A01", "A03", "A02")
  )

  drug_response_obj <- drug.response(response, lab_measurements, drug_purchases, c(-1, 0), c(0.1, 1))

  # Define an output prefix
  out_prefix <- "test_upset"

  # Run the function
  summarize_drug_purchases_upset(drug_response_obj, out_prefix)

  # Check that the output file was created
  output_file <- paste0(out_prefix, "_upset_plot.pdf")
  expect_true(file.exists(output_file))

  # Clean up the created file
  if (file.exists(output_file)) {
    file.remove(output_file)
  }
})

test_that("get_drug_purchases handles input correctly", {
  conn <- fg_data_connection(list(pheno = data.frame(
    FINNGENID = c("FG1", "FG2", "FG3"),
    SOURCE = c("PURCH", "PURCH", "PURCH"),
    CODE1 = c("A01", "A02", "A02"),
    CODE2 = c("", "", ""),
    CODE3 = c("", "", ""),
    CODE4 = c("1", "1", "1"),
    EVENT_AGE = c(21.0, 20.0, 35)
  )))

  expect_error(get_drug_purchases(conn, "A01", finngen_ids = "FG1"), NA)
})

# Helper function to create a realistic drug.response object for testing
create_test_drug_response_object <- function() {
  # Mock lab measurements
  mock_labs <- data.frame(
    FINNGENID = rep(paste0("FG", 1:4), each = 10),
    OMOP_CONCEPT_ID = 12345,
    EVENT_AGE = rep(c(30.1, 30.2, 30.8, 30.9, 31.0, 31.1, 31.2, 31.8, 31.9, 32.0), 4),
    MEASUREMENT_VALUE_HARMONIZED = rnorm(40, mean = 100, sd = 15)
  )

  # Mock drug purchases
  mock_purchases <- data.frame(
    FINNGENID = paste0("FG", 1:4),
    ATC = rep(c("DRUG_A", "DRUG_B"), each = 2),
    EVENT_AGE = c(31.05, 31.05, 31.05, 31.05),
    CODE1 = rep(c("DRUG_A", "DRUG_B"), each = 2) # for get_drug_purchases
  )

  # Add a clear outlier for testing the removal logic
  mock_labs$MEASUREMENT_VALUE_HARMONIZED[5] <- 500

  # Mimic the data preparation steps from create_drug_response
  dr_first_purchase <- mock_purchases %>%
    group_by(FINNGENID) %>%
    summarise(first_drug_age = first(EVENT_AGE), first_drug = first(ATC))

  all_measurements <- left_join(mock_labs, dr_first_purchase, by = "FINNGENID") %>%
    mutate(time_to_drug = .data$first_drug_age - .data$EVENT_AGE)

  # Create a minimal but valid drug.response object
  drug.response(
    responses = data.frame(FINNGENID = paste0("FG", 1:4)), # dummy data
    lab_measurements = all_measurements,
    drug_purchases = mock_purchases,
    before_period = c(0.01, 1.0), # years (time_to_drug is positive)
    after_period = c(-1.0, -0.01)  # years (time_to_drug is negative)
  )
}

test_that("plot_lab_value_distribution handles input and options correctly", {
  drug_resp_obj <- create_test_drug_response_object()

  # Test 1: Function returns a ggplot object without outlier removal
  p1 <- plot_lab_value_distribution(drug_resp_obj, remove_outliers = FALSE)
  expect_true(ggplot2::is_ggplot(p1))

  # Check that the outlier is present in the plot data
  expect_true(500 %in% p1$data$MEASUREMENT_VALUE_HARMONIZED)

  # Test 2: Function returns a ggplot object with outlier removal
  p2 <- plot_lab_value_distribution(drug_resp_obj, remove_outliers = TRUE)
  expect_true(ggplot2::is_ggplot(p2))

  # Check that the outlier has been removed from the plot data
  expect_false(500 %in% p2$data$MEASUREMENT_VALUE_HARMONIZED)

  # Test 3: Plot data has the correct structure
  p1_build <- ggplot2::ggplot_build(p1)
  expect_true("Before" %in% p1$data$period)
  expect_true("After" %in% p1$data$period)

  # Test 4: Function throws an error for invalid input
  expect_error(
    plot_lab_value_distribution(list()),
    "Input must be a drug.reponse object."
  )
})

test_that("create_drug_response handles covariates correctly", {
  # Mock minimal data
  mock_labs <- data.frame(FINNGENID = "FG1", EVENT_AGE = 30, MEASUREMENT_VALUE_HARMONIZED = 100, OMOP_CONCEPT_ID = "L1")
  mock_phenos <- data.frame(FINNGENID = "FG1", EVENT_AGE = 31, CODE1 = "D1", SOURCE = "PURCH")
  mock_covariates <- data.frame(
    FINNGENID = "FG1",
    SEX = "female", # Test "female" string
    AGE_AT_DEATH = 80,
    UNRELATED_COL = "data"
  )

  # Test 1: Add covariates successfully
  response_obj <- create_drug_response(
    kanta = mock_labs,
    phenos = mock_phenos,
    lablist = "L1",
    druglist = "D1",
    before_period = c(0.5, 1.5),
    after_period = c(-0.5, 0.5),
    covariates = mock_covariates,
    covariate_cols = c("SEX", "AGE_AT_DEATH")
  )

  # Check that columns were added to both dataframes in the object
  expect_true(all(c("SEX", "AGE_AT_DEATH") %in% colnames(response_obj$responses)))
  expect_true(all(c("SEX", "AGE_AT_DEATH") %in% colnames(response_obj$all_measurements)))

  # Check that unrequested columns were not added
  expect_false("UNRELATED_COL" %in% colnames(response_obj$responses))

  # Test 2: Function runs normally when covariates are not provided
  response_obj_no_cov <- create_drug_response(
    kanta = mock_labs,
    phenos = mock_phenos,
    lablist = "L1",
    druglist = "D1",
    before_period = c(0.5, 1.5),
    after_period = c(-0.5, 0.5)
  )
  expect_false("SEX" %in% colnames(response_obj_no_cov$responses))

  # Test 3: Function throws an error for missing covariate columns
  expect_error(
    create_drug_response(
      kanta = mock_labs,
      phenos = mock_phenos,
      lablist = "L1",
      druglist = "D1",
      before_period = c(0.5, 1.5),
      after_period = c(-0.5, 0.5),
      covariates = mock_covariates,
      covariate_cols = c("SEX", "NON_EXISTENT_COL") # One valid, one invalid
    ),
    "The following `covariate_cols` are not in the `covariates` dataframe: NON_EXISTENT_COL"
  )
})

# Helper function to create test data for BLUP analysis
create_test_blup_data <- function() {
  set.seed(123)  # For reproducibility

  # Create longitudinal lab measurements for 20 individuals
  n_individuals <- 20
  individuals <- paste0("FG", 1:n_individuals)

  # Generate multiple measurements per individual
  lab_data <- do.call(rbind, lapply(1:n_individuals, function(i) {
    n_measurements <- sample(3:8, 1)  # Random number of measurements
    base_age <- runif(1, 30, 50)
    ages <- base_age + sort(runif(n_measurements, 0, 10))

    # Individual-specific intercept and slope
    intercept <- rnorm(1, 100, 10)
    slope <- rnorm(1, -0.5, 0.2)  # Negative slope (decline over age)

    data.frame(
      FINNGENID = rep(individuals[i], n_measurements),
      OMOP_CONCEPT_ID = rep(c("3004410", "3023602")[((i-1) %% 2) + 1], n_measurements),
      EVENT_AGE = ages,
      MEASUREMENT_VALUE_HARMONIZED = intercept + slope * ages + rnorm(n_measurements, 0, 5),
      n = 1,
      first_drug_age = rep(base_age + 5, n_measurements),
      first_drug = rep("A10BJ", n_measurements),
      time_to_drug = rep(base_age + 5, n_measurements) - ages
    )
  }))

  # Create drug purchases data
  drug_purchases <- data.frame(
    FINNGENID = individuals,
    ATC = rep("A10BJ", n_individuals),
    EVENT_AGE = runif(n_individuals, 35, 55)
  )

  # Create sex data
  sex_data <- data.frame(
    FINNGENID = individuals,
    SEX = sample(c("male", "female"), n_individuals, replace = TRUE)
  )

  # Create drug.response object
  drug.response(
    responses = data.frame(FINNGENID = individuals),
    lab_measurements = lab_data,
    drug_purchases = drug_purchases,
    before_period = c(-1, 0),
    after_period = c(1/12, 1)
  )
}

test_that("calculate_blup_slopes works correctly", {
  # Create test data
  drug_resp_obj <- create_test_blup_data()
  sex_data <- data.frame(
    FINNGENID = paste0("FG", 1:20),
    SEX = sample(c("male", "female"), 20, replace = TRUE)
  )

  # Create temporary directory for output
  temp_dir <- tempdir()

  # Test 1: Function runs without error and produces output files
  blup_results <- calculate_blup_slopes(drug_resp_obj, sex_data = sex_data,
                                        output_dir = temp_dir)

  expect_true(is.list(blup_results))
  expect_true(length(blup_results) > 0)

  # Check that output files were created
  for (concept_id in names(blup_results)) {
    output_file <- file.path(temp_dir, paste0(concept_id, "_DF13.tsv"))
    expect_true(file.exists(output_file))

    # Read and check the output file
    output_data <- read.table(output_file, header = TRUE, sep = "\t")
    expect_equal(colnames(output_data), c("FID", "IID", paste0(concept_id, "_slope")))
    expect_equal(output_data$FID, output_data$IID)
  }

  # Test 2: Function works without sex data
  blup_results_no_sex <- calculate_blup_slopes(drug_resp_obj, sex_data = NULL,
                                                output_dir = temp_dir)
  expect_true(is.list(blup_results_no_sex))

  # Test 3: Function handles minimum measurements requirement
  blup_results_strict <- calculate_blup_slopes(drug_resp_obj, sex_data = sex_data,
                                                output_dir = temp_dir,
                                                min_measurements = 5)
  # Should have fewer individuals in the analysis
  for (concept_id in names(blup_results_strict)) {
    expect_true(blup_results_strict[[concept_id]]$n_individuals <=
                blup_results[[concept_id]]$n_individuals)
  }

  # Test 4: Function throws error for invalid input
  expect_error(
    calculate_blup_slopes(list()),
    "Input must be a drug.reponse object."
  )

  # Clean up temporary files
  for (concept_id in c("3004410", "3023602")) {
    file.remove(file.path(temp_dir, paste0(concept_id, "_DF13.tsv")))
  }
})

test_that("summarize_blup_results works correctly", {
  # Create test data and run BLUP analysis
  drug_resp_obj <- create_test_blup_data(include_sex = TRUE)
  temp_dir <- tempdir()

  blup_results <- calculate_blup_slopes(drug_resp_obj, output_dir = temp_dir)

  # Test summarize function
  summary_df <- summarize_blup_results(blup_results)

  expect_true(is.data.frame(summary_df))
  expect_equal(colnames(summary_df),
               c("OMOP_CONCEPT_ID", "n_individuals", "mean_slope",
                 "sd_slope", "min_slope", "max_slope"))
  expect_equal(nrow(summary_df), length(blup_results))

  # Check that slopes are negative on average (decline over age)
  expect_true(all(summary_df$mean_slope < 0))

  # Test with empty results
  empty_summary <- summarize_blup_results(list())
  expect_equal(nrow(empty_summary), 0)

  # Clean up
  for (concept_id in names(blup_results)) {
    file.remove(file.path(temp_dir, paste0(concept_id, "_DF13.tsv")))
  }
})